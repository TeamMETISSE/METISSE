BASE_dir=..

METISSE_src_dir=$(BASE_dir)/src

#compiler
FC=gfortran

WARN= -Wno-tabs -Wno-maybe-uninitialized -fimplicit-none  -Wsurprising  -Waliasing -pedantic
DEBUG= -g3 -O2 -fbacktrace -fcheck=all -fbounds-check

FFLAGS1= -ffree-form -fimplicit-none $(WARN) $(DEBUG)
FFLAGS2= -ffixed-form -fimplicit-none $(DEBUG)

# metisse in standalone mode
METISSE=metisse

SRCE0 = \
main_metisse.f90 evolv_metisse.f90 \
assign_commons_main.f90 assign_commons_other.f90 \
track_support.f90 z_support.f90 sse_support.f90 \
remnant_support.f90 interp_support.f90 \
METISSE_miscellaneous.f90 METISSE_zcnsts.f90 METISSE_star.f90 \
METISSE_hrdiag.f90 METISSE_deltat.f90

SRCE00 = \
assign_remnant.f hrdiag_remnant.f mrenv.f \

	
OBJT0 = $(SRCE0:.f90=.o)
OBJT00 = $(SRCE00:.f=.o)
	
track_support.o: $(METISSE_src_dir)/track_support.f90
z_support.o: $(METISSE_src_dir)/z_support.f90 track_support.o
interp_support.o: $(METISSE_src_dir)/interp_support.f90 track_support.o
sse_support.o: $(METISSE_src_dir)/sse_support.f90 track_support.o
remnant_support.o: $(METISSE_src_dir)/remnant_support.f90 track_support.o sse_support.o

assign_commons.o: $(METISSE_src_dir)/assign_commons_bse.f90 track_support.o remnant_support.o
METISSE_miscellaneous.o: $(METISSE_src_dir)/METISSE_miscellaneous.f90 track_support.o
comenv_lambda.o: $(METISSE_src_dir)/comenv_lambda.f90 track_support.o

METISSE_mlwind.o: $(METISSE_src_dir)/METISSE_mlwind.f90 track_support.o
METISSE_deltat.o: $(METISSE_src_dir)/METISSE_deltat.f90 track_support.o
METISSE_zcnsts.o: $(METISSE_src_dir)/METISSE_zcnsts.f90 track_support.o z_support.o
METISSE_star.o: $(METISSE_src_dir)/METISSE_star.f90 track_support.o interp_support.o sse_support.o
METISSE_hrdiag.o: $(METISSE_src_dir)/METISSE_hrdiag.f90 track_support.o interp_support.o sse_support.o remnant_support.o
METISSE_gntage.o: $(METISSE_src_dir)/METISSE_gntage.f90 track_support.o interp_support.o

evolv_metisse.o: $(METISSE_src_dir)/evolv_metisse.f90 track_support.o sse_support.o
main_metisse.o: $(METISSE_src_dir)/main_metisse.f90 track_support.o z_support.o
assign_commons_main.o: $(METISSE_src_dir)/assign_commons_main.f90 track_support.o remnant_support.o

%.o : $(METISSE_src_dir)/%.f90
	$(FC) $(FFLAGS1) -c $<
	
%.o : $(METISSE_src_dir)/dummies/%.f
	$(FC) $(FFLAGS2) -c $<


$(METISSE): $(OBJT0) $(OBJT00)
	$(FC) $(OMP) $(FFLAGS1) -o $(BASE_dir)/$(METISSE) $^

all: $(METISSE)


clean:
	rm -f $(BASE_dir)/$(METISSE)
	rm -f *.mod *.o

