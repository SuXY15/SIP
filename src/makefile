make:
	ifort -O3 SIP_Spherical.f90 -o SIP
clean:
	rm *.txt
flame:
	python3 ./tools/showData.py flame_position.txt
show:
	python3 ./tools/showData.py data${data}.txt
display:
	python3 ./tools/dispData.py ${jump}
save:
	mv *.txt ${dir}