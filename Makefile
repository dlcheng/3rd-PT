#--------------------------------------- Select target computer
#SYSTYPE="Workstation"
SYSTYPE="Mac"

#--------------------------------------- Adjust settings for target computer
ifeq ($(SYSTYPE),"Workstation")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall 
GSL_INCL =  -I/home/dlcheng/Install/gsl/include
GSL_LIBS =  -L/home/dlcheng/Install/gsl/lib
endif

ifeq ($(SYSTYPE),"Mac")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall
GSL_LIBS=   -L/usr/local/Cellar/gsl/1.16/lib  -Wl
GSL_INCL =  -I/usr/local/Cellar/gsl/1.16/include
endif

OPTIONS =  $(OPTIMIZE) $(OPT) 

EXEC   =  3rd_PT

OBJS   =  allvars.o growth.o main.o kernel_intg.o spectra.o set_params.o tk_bbks.o tk_from_file.o linear_p.o lp_from_file.o norm_p.o

INCL   = allvars.h proto.h define.h Makefile

CFLAGS = $(OPTIONS) $(GSL_INCL) 

LIBS   = $(GSL_LIBS) -lgsl -lgslcblas -lm 

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS) $(EXEC) *.gch
