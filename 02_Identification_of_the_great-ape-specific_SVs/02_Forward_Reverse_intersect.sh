# For great apes (V38 CCP GGO PAB)
for i in V38 CCP GGO PAB 
do
bedtools intersect -a ../F-variants/NLE-"$i"-SV.txt-DEL -b ../R-variants/NLE-"$i"-SV.txt-DEL -wa -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > NLE-"$i"-SV.txt-DEL
bedtools intersect -a ../F-variants/NLE-"$i"-SV.txt-INS-plus -b ../R-variants/NLE-"$i"-SV.txt-INS-plus -wa -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > NLE-"$i"-SV.txt-INS-plus
done 

# For outgroup (RM10)
cat ../F-variants/NLE-RM10-SV.txt-DEL ../R-variants/NLE-RM10-SV.txt-DEL |sort -k1,1 -s -V -k2n,2  > NLE-RM10-SV.txt-DEL
cat ../F-variants/NLE-RM10-SV.txt-INS-plus ../R-variants/NLE-RM10-SV.txt-INS-plus |sort -k1,1 -s -V -k2n,2 > NLE-RM10-SV.txt-INS-plus
