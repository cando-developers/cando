
cando-run:
	$(HOME)/Development/amber16/bin/sander -O -i 01_Min.in -o 01_Min_cando.out -p alanine_cando.top -c alanine_cando.crd

leap-run:
	$(HOME)/Development/amber16/bin/sander -O -i 01_Min.in -o 01_Min_leap.out -p alanine_leap.top -c alanine_leap.crd

leap-run-no-box:
	$(HOME)/Development/amber16/bin/sander -O -i 01_Min.in -o 01_Min_leap_no_box.out -p alanine_leap_no_box.top -c alanine_leap_no_box.crd

cando-chemdraw-run:
	$(HOME)/Development/amber16/bin/sander -O -i 01_Min.in -o 01_Min_chemagg.out -p chemagg.top -c chemagg.crd -r chemagg.rst -inf chemagg.mdinfo

cando-water-run:
	$(HOME)/Development/amber16/bin/sander -O -i 01_Min.in -o 01_Min_chemagg_water.out -p chemagg_water.top -c chemagg_water.crd -r chemagg_water.rst -inf chemagg_water.mdinfo

cando-run-heat:
	$(HOME)/Development/amber16/bin/sander -O -i 02_Heat.in -o 02_Heat.out -p chemagg_water.top -c chemagg_water.rst -r 02_Heat.rst -x 02_Heat.mdcrd -inf 02_Heat.mdinfo


leap-water-run:
	$(HOME)/Development/amber16/bin/sander -O -i 01_Min_chemagg_water.in -o 01_Min_chemagg_water.out -p prmtop -c inpcrd -r leap_water.rst -inf leap_water.mdinfo

leap-run-heat:
	$(HOME)/Development/amber16/bin/sander -O -i 02_Heat.in -o 02_Heat.out -p prmtop -c leap_water.rst -r 02_Heat_leap.rst -x 02_Heat_leap.mdcrd -inf 02_Heat_leap.mdinfo
