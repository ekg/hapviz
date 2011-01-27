#!/usr/bin/perl

$color_A=32;
$color_C=34;
$color_G=33;
$color_T=31;

while(<>) {
	s/(A)/\033\[1;${color_A}m$1\033\[0m/ig;
	s/(C)/\033\[1;${color_C}m$1\033\[0m/ig;
	s/(G)/\033\[1;${color_G}m$1\033\[0m/ig;
	s/(T)/\033\[1;${color_T}m$1\033\[0m/ig;
	print;
}
