for file in ./*.int
do 
	echo $file
	grep "Interatomic Surface   " $file|cut -d: -f2 >> IS.txt
	#echo "                   " >> IS.txt
	grep "Flux of J_\*X(r)  " $file |cut -d\) -f2 |cut -du -f2 |cut -d= -f2|cut -dn -f1 >> JX.txt
	cp JX.txt foo.txt.tmp
	sed '$ d' foo.txt.tmp > JX.txt
	rm -f foo.txt.tmp
	
	grep "Flux of J_\*Y(r)  " $file |cut -d\) -f2 |cut -du -f2 |cut -d= -f2|cut -dn -f1 >> JY.txt
        cp JY.txt foo.txt.tmp
        sed '$ d' foo.txt.tmp > JY.txt
        rm -f foo.txt.tmp
	
	
	grep "Flux of J_\*Z(r)  " $file |cut -d\) -f2 |cut -du -f2 |cut -d= -f2|cut -dn -f1 >> JZ.txt
        cp JZ.txt foo.txt.tmp
        sed '$ d' foo.txt.tmp > JZ.txt
        rm -f foo.txt.tmp
	
	
	
	pr -mts' ' IS.txt JX.txt JY.txt JZ.txt >>TotalJ.txt	
	rm IS.txt JX.txt JY.txt JZ.txt

done



