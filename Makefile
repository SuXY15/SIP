PROGRAM := SIP
FC      := gfortran
FCFLAGS := -c -ffree-line-length-none
FCFLAGS += -g -fbacktrace -fno-align-commons -fbounds-check
FDFLAGS := -Ddebug
FLFLAGS := 

BINDIR := bin
SRCDIR := src
SRCS   := $(wildcard $(addsuffix /*.f90, $(SRCDIR)))
OBJS   := $(SRCS:$(SRCDIR)/%.f90=$(BINDIR)/%.o)
PRG_OBJ := $(BINDIR)/$(PROGRAM).o

# For compiling
default:
	mkdir -p bin
	mkdir -p data
	make clean 
	make $(PROGRAM)

$(OBJS): $(BINDIR)/%.o:$(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) -cpp $(FDFLAGS) -o $@ $<

$(PROGRAM):$(OBJS)
	$(FC) $(FLFLAGS) -o $@ $^
	mv *.mod $(BINDIR)

# Dependencies
$(BINDIR)/init.o: $(BINDIR)/para.o
$(BINDIR)/main.o: $(BINDIR)/init.o

supernova:
	export FDFLAGS="$(FDFLAGS) -Dsupernova"
	make default

create_case:
	mkdir data/${case}
	mv data/*.txt data/${case}
clean_case:
	rm data/${case}/*.dat

get_flame:
	python3 tools/getFlame.py data/${case}
show_flame:
	python3 tools/showFlame.py data/${case}/flame_position.dat ${start} ${end}

show_data:
	python3 tools/showData.py data/${case}/${data}.dat
display_data:
	python3 tools/dispData.py data/${case} ${jump} 

test:
	python3 tools/test.py
clean:
	-rm -f $(OBJS) $(PROGRAM) *.mod

