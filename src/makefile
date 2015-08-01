

SRC-SUB-DIRS = geom units adapt chem # kinematics

sub-prebuild:
	$(foreach var,$(SRC-SUB-DIRS), make -C $(var)/ $@; )
