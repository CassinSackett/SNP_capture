BEGIN{
	FS="\t"
}
NR==1 {
	name=$x;
	print (name) > (name ".homs");
	next
}
{
split($x, a, "/"); if (a[1] == a[2]) print ($1 "\t" $2 "\t" $x) >> (name ".homs")
	}
END {
}
