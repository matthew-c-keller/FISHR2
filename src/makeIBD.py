f = open("parse_bmatch_Beagle.txt","r")
fw= open("InputIBD_Beagle.txt","w")
fr= f.readlines()
for each in fr:
	each = each.split("\t")
	p1 = each[0].split(" ")[0]
	p2 = each[1].split(" ")[0]
	bp1 = each[3].split(" ")[0]
	bp2 = each[3].split(" ")[1]
	output = "{} {} {} {}\n".format(p1,p2,bp1,bp2)
	fw.write(output)