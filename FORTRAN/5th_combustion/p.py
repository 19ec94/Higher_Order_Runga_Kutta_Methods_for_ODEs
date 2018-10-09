file=open("sol_exact.out","r")
fileout=open("sol_exact_withoutgap.out","w")
for line in file:
	print(line.lstrip("   ").rstrip("\n"),file=fileout)
