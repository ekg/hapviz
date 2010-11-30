#include <string>
#include <iostream>

using namespace std;

class IndelAllele {
    friend ostream& operator<<(ostream&, const IndelAllele&);
    friend bool operator==(const IndelAllele&, const IndelAllele&);
    friend bool operator!=(const IndelAllele&, const IndelAllele&);
    friend bool operator<(const IndelAllele&, const IndelAllele&);
public:
    bool insertion;
    int length;
    int position;
    string sequence;

    IndelAllele(bool i, int l, int p, string s)
        : insertion(i), length(l), position(p), sequence(s)
    { }
};

ostream& operator<<(ostream& out, IndelAllele& indel) {
    string t = indel.insertion ? "i" : "d";
    out << t <<  ":" << indel.position << ":" << indel.sequence;
    return out;
}

bool operator==(const IndelAllele& a, const IndelAllele& b) {
    return (a.insertion == b.insertion
            && a.length == b.length
            && a.position == b.position
            && a.sequence == b.sequence);
}

bool operator!=(const IndelAllele& a, const IndelAllele& b) {
    return !(a==b);
}

bool operator<(const IndelAllele& a, const IndelAllele& b) {
    if (a == b)
        return false;
    if (!a.insertion && b.insertion) {
        return true;
    } else if (a.position < b.position) {
        return true;
    } else if (a.length < b.length) {
        return true;
    } else if (a.length == b.length) { 
        return a.sequence < b.sequence;
    } else {
        return false;
    }
}
