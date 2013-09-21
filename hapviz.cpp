#include <iostream>
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <signal.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include "BamAlignment.h"
#include "BamMultiReader.h"
#include "BamReader.h"
#include "BamWriter.h"
#include "Fasta.h"
#include "SmithWatermanGotoh.h"
#include "BandedSmithWaterman.h"

#include "IndelAllele.h"

#include "levenshtein.h"

using namespace std;
using namespace BamTools;



// drawn from bamtools API
//
const char REVCOMP_LOOKUP[] = {'T', 0, 'G', 'H', 0, 0, 'C', 'D', 0, 0, 0, 0, 'K', 'N', 0, 0, 0, 'Y', 'W', 'A', 'A', 'B', 'S', 'X', 'R', 0 };

void ReverseComplement(std::string& sequence) {
    
    // do complement
    size_t seqLength = sequence.length();
    for ( size_t i = 0; i < seqLength; ++i )
        sequence.replace(i, 1, 1, REVCOMP_LOOKUP[(int)sequence.at(i) - 65]);
    
    // reverse it
    reverse(sequence.begin(), sequence.end());
}

// print BamAlignment in FASTQ format
// N.B. - uses QueryBases NOT AlignedBases
void writeFastq(ostream& out, const BamAlignment& a) { 
 
    // @BamAlignment.Name
    // BamAlignment.QueryBases
    // +
    // BamAlignment.Qualities
    //
    // N.B. - QueryBases are reverse-complemented (& Qualities reversed) if aligned to reverse strand .
    //        Name is appended "/1" or "/2" if paired-end, to reflect which mate this entry is.
  
    // handle paired-end alignments
    string name = a.Name;
    if ( a.IsPaired() )
        name.append( (a.IsFirstMate() ? "/1" : "/2") );
  
    // handle reverse strand alignment - bases & qualities
    string qualities = a.Qualities;
    string sequence  = a.QueryBases;
    if ( a.IsReverseStrand() ) {
        ReverseComplement(sequence);
    }
  
    // write to output stream
    out << "@" << name << endl
        << sequence    << endl
        << "+"         << endl
        << qualities   << endl;
}

short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}

bool sufficientQuality(BamAlignment& alignment, int minbaseq) {
    for (string::iterator q = alignment.Qualities.begin(); q != alignment.Qualities.end(); ++q) {
        if (!(qualityChar2ShortInt(*q) >= minbaseq))
            return false;
    }
    return true;
}

string indelClassStr(int indelClass) {
    stringstream s;
    if (indelClass == 0)
        return "M";
    if (indelClass > 0)
        s << "I" << abs(indelClass);
    if (indelClass < 0)
        s << "D" << abs(indelClass);
    return s.str();
}

// directly pulled from bamtools source
// Parses a region string, does validation (valid ID's, positions), stores in Region struct
// Returns success (true/false)
bool ParseRegionString(const std::string& regionString, const BamMultiReader& reader, BamRegion& region) {
  
    // -------------------------------
    // parse region string
  
    // check first for empty string
    if ( regionString.empty() ) 
        return false;   
    
    // non-empty string, look for a colom
    size_t foundFirstColon = regionString.find(':');
    
    // store chrom strings, and numeric positions
    string startChrom;
    string stopChrom;
    int startPos;
    int stopPos;
    
    // no colon found
    // going to use entire contents of requested chromosome 
    // just store entire region string as startChrom name
    // use BamReader methods to check if its valid for current BAM file
    if ( foundFirstColon == string::npos ) {
        startChrom = regionString;
        startPos   = 0;
        stopChrom  = regionString;
        stopPos    = -1;
    }
    
    // colon found, so we at least have some sort of startPos requested
    else {
      
        // store start chrom from beginning to first colon
        startChrom = regionString.substr(0,foundFirstColon);
        
        // look for ".." after the colon
        size_t foundRangeDots = regionString.find("..", foundFirstColon+1);
        
        // no dots found
        // so we have a startPos but no range
        // store contents before colon as startChrom, after as startPos
        if ( foundRangeDots == string::npos ) {
            startPos   = atoi( regionString.substr(foundFirstColon+1).c_str() ); 
            stopChrom  = startChrom;
            stopPos    = -1;
        } 
        
        // ".." found, so we have some sort of range selected
        else {
          
            // store startPos between first colon and range dots ".."
            startPos = atoi( regionString.substr(foundFirstColon+1, foundRangeDots-foundFirstColon-1).c_str() );
          
            // look for second colon
            size_t foundSecondColon = regionString.find(':', foundRangeDots+1);
            
            // no second colon found
            // so we have a "standard" chrom:start..stop input format (on single chrom)
            if ( foundSecondColon == string::npos ) {
                stopChrom  = startChrom;
                stopPos    = atoi( regionString.substr(foundRangeDots+2).c_str() );
            }
            
            // second colon found
            // so we have a range requested across 2 chrom's
            else {
                stopChrom  = regionString.substr(foundRangeDots+2, foundSecondColon-(foundRangeDots+2));
                stopPos    = atoi( regionString.substr(foundSecondColon+1).c_str() );
            }
        }
    }

    // -------------------------------
    // validate reference IDs & genomic positions
    
    const RefVector references = reader.GetReferenceData();
    
    // if startRefID not found, return false
    int startRefID = reader.GetReferenceID(startChrom);
    if ( startRefID == (int)references.size() ) return false;  
    
    // if startPos is larger than reference, return false
    const RefData& startReference = references.at(startRefID);
    if ( startPos > startReference.RefLength ) return false;
    
    // if stopRefID not found, return false
    int stopRefID = reader.GetReferenceID(stopChrom);
    if ( stopRefID == (int)references.size() ) return false;
    
    // if stopPosition larger than reference, return false
    const RefData& stopReference = references.at(stopRefID);
    if ( stopPos > stopReference.RefLength ) return false;
    
    // if no stopPosition specified, set to reference end
    if ( stopPos == -1 ) stopPos = stopReference.RefLength;  
    
    // -------------------------------
    // set up Region struct & return
    
    region.LeftRefID = startRefID;
    region.LeftPosition = startPos;
    region.RightRefID = stopRefID;;
    region.RightPosition = stopPos;
    return true;
}

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] bam file [... bam file]" << endl
         << endl
         << "options:" << endl 
         << "    -b, --bam           bam file to from which to extract reads" << endl
         << "                        (may also be specified as positional arguments)" << endl
         << "    -r, --region        region from which to extract grouped reads" << endl
        //<< "    -v, --visualize     Print an ASCII-art style pileup for each read-haplotype group" << endl
         << "    -a, --show-all      When visualizing, show all alignments, not just variant ones." << endl
        // << "    -o, --output-pefix  prefix for output bam files with grouped reads" << endl
         << "    -f, --reference     FASTA reference against which alignments have been aligned" << endl
         << "    -q, --min-base-quality" << endl
         << "                        minimum base quality required for all bases in a read" << endl
         << "    -H, --haplotypes    Write out haplotypes observed in the target region window" << endl
         << endl
         << "Displays haplotype groups from the specified region across the BAM files provided as input." << endl
         << endl;
    exit(0);
}


int main (int argc, char** argv) {

    int c;

    vector<string> bam_files;
    string output_prefix = "bamgroups";
    string region_str = "";
    int minbaseq = 0;
    bool hasRef = false;
    FastaReference reference;
    bool visualize = true;
    bool printHaplotypeCounts = false;
    bool realign = false;
    bool showAllReads = false;
    bool useStdin = false;

    if (argc == 1)
        printSummary(argv);

    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"visualize", no_argument, 0, 'v'},
                {"show-all", no_argument, 0, 'a'},
                {"bam",  required_argument, 0, 'b'},
                {"output-prefix",  required_argument, 0, 'o'},
                {"region", required_argument, 0, 'r'},
                {"min-base-quality", required_argument, no_argument, 'q'},
                {"reference", required_argument, 0, 'f'},
                {"haplotype-counts", no_argument, 0, 'H'},
                {"stdin", no_argument, 0, 'c'},
                {0, 0, 0, 0}
            };

        int option_index = 0;

        c = getopt_long (argc, argv, "hvHcab:f:o:r:q:",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c)
        {
        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf ("option %s", long_options[option_index].name);
            if (optarg)
                printf (" with arg %s", optarg);
            printf ("\n");
            break;

        case 'b':
            bam_files.push_back(optarg);
            break;
 
        case 'o':
            output_prefix = optarg;
            break;
 
        case 'r':
            region_str = optarg;
            break;

        case 'H':
            visualize = false;
            printHaplotypeCounts = true;
            break;

        case 'f':
            reference.open(optarg);
            hasRef = true;
            break;
 
        case 'q':
            minbaseq = atoi(optarg);
            break;

        case 'a':
            showAllReads = true;
            break;
 
        case 'v':
            visualize = true;
            break;

        case 'c':
            useStdin = true;
            bam_files.push_back("stdin");
            break;
 
        case 'h':
            printSummary(argv);
            exit(0);
            break;
          
        case '?':
            /* getopt_long already printed an error message. */
            printSummary(argv);
            exit(1);
            break;
 
        default:
            abort ();
        }
    }

    // use rest of command-line arguments as BAM files
    if (optind < argc) {
        while (optind < argc) {
            bam_files.push_back(argv[optind++]);
        }
    }


    if (!hasRef) {
        cerr << "you must supply a reference sequence (-f)" << endl;
        return 1;
    }

    if (bam_files.empty()) {
        useStdin = true;
        bam_files.push_back("stdin");
    }

    // open the listed BAM files
    BamMultiReader reader;
    if (!region_str.empty()) {
        if (!reader.Open(bam_files, true)) {
            cerr << "Could not open bam file(s)" << endl;
            exit(1);
        }
    } else {
        if (!reader.Open(bam_files, false)) {
            cerr << "Could not open bam file(s)" << endl;
            exit(1);
        }
    }

    // set target region if specified
    BamRegion region;
    if (!region_str.empty()) {
        if (region_str.find("-") != string::npos) {
            region_str.replace(region_str.find("-"), 1, "..");
        }
        ParseRegionString(region_str, reader, region);
        if (!reader.SetRegion(region)) {
            cerr << "Could not set region to " << region_str << endl;
            exit(1);
        }
    }

    map<int, string> referenceIDToName;
    // store the names of all the reference sequences in the BAM file
    vector<RefData> referenceSequences = reader.GetReferenceData();
    int i = 0;
    for (RefVector::iterator r = referenceSequences.begin(); r != referenceSequences.end(); ++r) {
        referenceIDToName[i] = r->RefName;
        ++i;
    }

    // retrieve header information
    map<string, string> readGroupToSampleNames;

    string bamHeader = reader.GetHeaderText();

    vector<string> headerLines = split(bamHeader, '\n');

    for (vector<string>::const_iterator it = headerLines.begin(); it != headerLines.end(); ++it) {

        // get next line from header, skip if empty
        string headerLine = *it;
        if ( headerLine.empty() ) { continue; }

        // lines of the header look like:
        // "@RG     ID:-    SM:NA11832      CN:BCM  PL:454"
        //                     ^^^^^^^\ is our sample name
        if ( headerLine.find("@RG") == 0 ) {
            vector<string> readGroupParts = split(headerLine, "\t ");
            string name = "";
            string readGroupID = "";
            for (vector<string>::const_iterator r = readGroupParts.begin(); r != readGroupParts.end(); ++r) {
                vector<string> nameParts = split(*r, ":");
                if (nameParts.at(0) == "SM") {
                    name = nameParts.at(1);
                } else if (nameParts.at(0) == "ID") {
                    readGroupID = nameParts.at(1);
                }
            }
            if (name == "") {
                cerr << " could not find SM: in @RG tag " << endl << headerLine << endl;
                return 1;
            }
            if (readGroupID == "") {
                cerr << " could not find ID: in @RG tag " << endl << headerLine << endl;
                return 1;
            }
            //string name = nameParts.back();
            //mergedHeader.append(1, '\n');
            //cerr << "found read group id " << readGroupID << " containing sample " << name << endl;
            readGroupToSampleNames[readGroupID] = name;
        }
    }

    // groups reads by the +/- variance they represent relative to the reference
    map<int, vector<BamAlignment> > indelGroupsByVariance;
    //  pair<ins   , del>
    typedef map<vector<IndelAllele>, vector<BamAlignment> >
        AlignmentAlleleGrouping;
    AlignmentAlleleGrouping indelGroupsByInsDelSeq;

    BamAlignment alignment;

    // extract the groups of alignments
    while (reader.GetNextAlignment(alignment)) {

        string referenceSequence = reference.getSubSequence(referenceIDToName[alignment.RefID], alignment.Position, alignment.GetEndPosition());

        // skip this alignment if we are not using duplicate reads (we remove them by default)
        // skip unmapped alignments, as they cannot be used in the algorithm
        // skip alignments which are non-primary
        if (alignment.IsDuplicate()
            || !alignment.IsMapped()
            || !alignment.IsPrimaryAlignment()) {
            continue;
        }

        // or... if we don't have sufficient quality in the alignment
        if (minbaseq > 0 && !sufficientQuality(alignment, minbaseq))
            continue;

        int rp = 0;  // read position, 0-based relative to read
        int csp = 0; // current sequence position, 0-based relative to currentSequence
        int sp = alignment.Position;  // sequence position

        vector<IndelAllele> indels;

        vector<CigarOp>::const_iterator cigarIter = alignment.CigarData.begin();
        vector<CigarOp>::const_iterator cigarEnd  = alignment.CigarData.end();
        for ( ; cigarIter != cigarEnd; ++cigarIter ) {
            unsigned int l = cigarIter->Length;
            char t = cigarIter->Type;

            if (t == 'M') { // match or mismatch

                sp += l;
                csp += l;
                rp += l;

            } else if (t == 'D') { // deletion

                indels.push_back(IndelAllele(false, l, sp, referenceSequence.substr(csp, l)));
                sp += l;  // update sample position
                csp += l;

            } else if (t == 'I') { // insertion

                indels.push_back(IndelAllele(true, l, sp, alignment.QueryBases.substr(rp, l)));
                rp += l;

                // handle other cigar element types
            } else if (t == 'S') { // soft clip, clipped sequence present in the read not matching the reference
                // skip these bases in the read
                /*
                  if (rp == 0) {
                  alignment.QueryBases = alignment.QueryBases.substr(l);
                  } else {
                  alignment.QueryBases = alignment.QueryBases.substr(0, alignment.QueryBases.size() - l);
                  }*/
                rp += l;// sp += l; csp += l;
            } else if (t == 'H') { // hard clip on the read, clipped sequence is not present in the read
            } else if (t == 'N') { // skipped region in the reference not present in read, aka splice
                sp += l; csp += l;
            }
        } // end cigar iter loop

        indelGroupsByInsDelSeq[indels].push_back(alignment);


        /*
          int variance = 0;

          for (vector<CigarOp>::const_iterator c = alignment.CigarData.begin(); c != alignment.CigarData.end(); ++c) {
          if (c->Type == 'D')
          variance -= c->Length;
          if (c->Type == 'I')
          variance += c->Length;
          }
          indelGroupsByVariance[variance].push_back(alignment);
        */
    }


    // now the alignments are grouped by reference-relative insertion and
    // deletion sequence within the window
    // the next step is to find the biggest group(s)
    // establish the haplotypes implied by these groups
    // realign all the reads against the haplotypes
    // and take the best (least-variant) alignment for each read
    
    if (visualize || printHaplotypeCounts) {
        for (AlignmentAlleleGrouping::iterator g = indelGroupsByInsDelSeq.begin();
             g != indelGroupsByInsDelSeq.end(); ++g) {

            vector<BamAlignment>& alignments = g->second;
            vector<IndelAllele> indels = g->first;

            stringstream indelsstream;
            string indelseq;
            bool inRegion = false;

            // ignore reads without indels
            if (!showAllReads && indels.empty()) {
                continue;
            } else {
                inRegion = true;
            }

            for (vector<IndelAllele>::iterator i = indels.begin(); i != indels.end(); ++i) {
                indelsstream << "," << *i;
                if (!region_str.empty()) {
                    if (i->position >= region.LeftPosition && i->position < region.RightPosition) {
                        inRegion = true;
                    }
                }
            }
            if (region_str.empty())
                inRegion = true;
            if (!inRegion)
                continue;

            // indelseq is the indel allele specifier
            indelseq = indelsstream.str();
            if (!indelseq.empty())
                indelseq = indelseq.substr(1);

            cout << alignments.size() << " " << indelseq << endl;
            BamAlignment& firstAlignment = alignments.front();
            int firstpos = firstAlignment.Position;
            int lastpos = firstAlignment.GetEndPosition();
            for (vector<BamAlignment>::iterator a = alignments.begin(); a != alignments.end(); ++a) {
                if (a->GetEndPosition() > lastpos)
                    lastpos = a->GetEndPosition();
            }
            int offset = 0;
            string refseq = reference.getSubSequence(referenceIDToName[firstAlignment.RefID], firstpos, lastpos - firstpos + 1);
            for (vector<IndelAllele>::iterator i = indels.begin(); i != indels.end(); ++i) {
                if (i->insertion) {
                    refseq.insert(i->position - firstAlignment.Position + offset, string(i->length, '-'));
                    offset += i->length;
                }
            }

            unsigned int maxReadNameSize = 0;
            for (vector<BamAlignment>::iterator a = alignments.begin(); a != alignments.end(); ++a) {
                if (a->Name.size() > maxReadNameSize) {
                    maxReadNameSize = a->Name.size();
                }
            }

            stringstream refposs;
            refposs << firstpos;
            string refpos = refposs.str();
            if (visualize) {
                cout << string(5 + maxReadNameSize - refpos.size(), ' ') << refpos << "   " << refseq << "   " << lastpos << endl;
            }

            map<string, vector<string> > alignmentsBySample;
            for (vector<BamAlignment>::iterator a = alignments.begin(); a != alignments.end(); ++a) {
                stringstream cigar;
                for (vector<CigarOp>::const_iterator cigarIter = a->CigarData.begin();
                     cigarIter != a->CigarData.end(); ++cigarIter) {
                    cigar << cigarIter->Length << cigarIter->Type;
                }
                //string cigarstr = cigar.str();
                string readGroup;
                a->GetTag("RG", readGroup);
                int pos = a->Position;
                string& samplename = readGroupToSampleNames[readGroup];
                // hacky...
                int pad = (maxReadNameSize) - a->Name.size() +4;
                if (pad < 0) pad = 1;
                stringstream alstr;
                alstr << a->Name << (a->IsFirstMate() ? ".1" : ".2") << (a->IsReverseStrand() ? " -" : " +") << string(pad, ' ')
                      << string(pos - firstpos, ' ') << a->AlignedBases;
                alignmentsBySample[samplename].push_back(alstr.str());
            }

            if (visualize) {
                for (map<string, vector<string> >::iterator s = alignmentsBySample.begin(); s != alignmentsBySample.end(); ++s) {
                    const string& sampleName = s->first;
                    vector<string>& aligns = s->second;
                    cout << endl;
                    cout << sampleName << endl;
                    for (vector<string>::iterator a = aligns.begin(); a != aligns.end(); ++a) {
                        cout << *a << endl;
                    }
                }
            }

            cout << endl;

        }
    }

    if (realign) {

        AlignmentAlleleGrouping::iterator biggestGroup = indelGroupsByInsDelSeq.begin();

        for (AlignmentAlleleGrouping::iterator g = indelGroupsByInsDelSeq.begin();
             g != indelGroupsByInsDelSeq.end(); ++g) {

            vector<BamAlignment>& alignments = g->second;
            vector<IndelAllele> indels = g->first;

            // ignore reads without indels
            if (indels.empty()) {
                continue;
            } else if (biggestGroup == indelGroupsByInsDelSeq.begin()) {
                biggestGroup = g;
            }
            stringstream indelsstream;
            string indelseq;
            bool inRegion = false;
            for (vector<IndelAllele>::iterator i = indels.begin(); i != indels.end(); ++i) {
                indelsstream << "," << *i;
                if (!region_str.empty()) {
                    if (i->position >= region.LeftPosition && i->position < region.RightPosition) {
                        inRegion = true;
                    }
                }
            }
            if (region_str.empty())
                inRegion = true;
            if (!inRegion)
                continue;

            indelseq = indelsstream.str();
            if (!indelseq.empty())
                indelseq = indelseq.substr(1);

            if (alignments.size() > biggestGroup->second.size()) {
                biggestGroup = g;
            }

        }

        cerr << "biggest group size " << biggestGroup->second.size() << endl;
        // establish the haplotype of this group

        vector<BamAlignment>& alignments = biggestGroup->second;
        vector<IndelAllele> indels = biggestGroup->first;

        BamAlignment& firstAlignment = alignments.front();
        int firstpos = firstAlignment.Position;
        int lastpos = firstAlignment.GetEndPosition();
        for (vector<BamAlignment>::iterator a = alignments.begin(); a != alignments.end(); ++a) {
            if (a->GetEndPosition() > lastpos)
                lastpos = a->GetEndPosition();
        }
        int offset = 0;
        string haplotype = reference.getSubSequence(referenceIDToName[firstAlignment.RefID], firstpos, lastpos - firstpos);
        for (vector<IndelAllele>::iterator i = indels.begin(); i != indels.end(); ++i) {
            if (i->insertion) {
                haplotype.insert(i->position - firstAlignment.Position + offset, i->sequence);
                offset += i->length;
            } else {
                haplotype.erase(i->position - firstAlignment.Position + offset, i->length);
                offset -= i->length;
            }
        }

        cerr << haplotype << endl;
    }

    /*

	// initialize
	
	const unsigned int referenceLen = strlen(reference);
	const unsigned int queryLen     = strlen(query);

	unsigned int referencePos;
	string cigar;

	// create a new Smith-Waterman alignment object
    if (bandwidth > 0) {
    pair< pair<unsigned int, unsigned int>, pair<unsigned int, unsigned int> > hr;
    hr.first.first   = 5;
    hr.first.second  = 13;
    hr.second.first  = 0;
    hr.second.second = 8;
    CBandedSmithWaterman bsw(matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty, bandwidth);
    bsw.Align(referencePos, cigar, reference, referenceLen, query, queryLen, hr);
    } else {
    CSmithWatermanGotoh sw(matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
    sw.Align(referencePos, cigar, reference, referenceLen, query, queryLen);
    }

    printf("%s %3u\n", cigar.c_str(), referencePos);

    */

    /*
      if (output_fastq) {

      for (map<int, vector<BamAlignment> >::iterator va = indelGroups.begin(); va != indelGroups.end(); ++va) {
      int indelClass = va->first;
      vector<BamAlignment>& alignments = va->second;
      ofstream outputFile;
      outputFile.open((output_prefix + "." + indelClassStr(indelClass) + ".fq").c_str());
      for (vector<BamAlignment>::iterator a = alignments.begin(); a != alignments.end(); ++a) {
      writeFastq(outputFile, *a);
      }
      outputFile.close();
      }

      } else {

      const string headerText = reader.GetHeaderText();
      RefVector references = reader.GetReferenceData();

      for (map<int, vector<BamAlignment> >::iterator va = indelGroups.begin(); va != indelGroups.end(); ++va) {
      int indelClass = va->first;
      vector<BamAlignment>& alignments = va->second;
      BamWriter writer;
      writer.Open(output_prefix + "." + indelClassStr(indelClass) + ".bam", headerText, references);
      for (vector<BamAlignment>::iterator a = alignments.begin(); a != alignments.end(); ++a) {
      writer.SaveAlignment(*a);
      }
      writer.Close();
      }

      }
    */

    reader.Close();

    return 0;

}
