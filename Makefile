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

clean:
	rm -f $(OBJS) $(PROGRAM) *.mod

visual_data:
	python3 tools/visual_data.py data/${case}/${data}.dat
visual_flame:
	python3 tools/visual_flame.py data/${case}/flame_position.dat ${start} ${end}
visual_evolve:
	python3 tools/visual_evolve.py data/${case} ${jump} 
visual_getFlame:
	python3 tools/visual_getFlame.py data/${case}

