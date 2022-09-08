#-- for now these MUST point to the included "samtools-0.x.x" and "gclib" sub-directories
HTSLIB  := ./htslib
#-- 
LIBDEFLATE := ${HTSLIB}/xlibs/lib/libdeflate.a
LIBBZ2 := ${HTSLIB}/xlibs/lib/libbz2.a
LIBLZMA := ${HTSLIB}/xlibs/lib/liblzma.a

GDIR := ./gclib
#--
INCDIRS := -I. -I${GDIR} -I${HTSLIB}

CXX   := $(if $(CXX),$(CXX),g++)

# BASEFLAGS := -Wall -Wextra ${INCDIRS} ${PYTHONFLAGS} -fsigned-char -D_FILE_OFFSET_BITS=64 \
# -D_LARGEFILE_SOURCE -std=c++11 -fno-strict-aliasing -fno-rtti

BASEFLAGS := -Wall -Wextra ${INCDIRS} -fsigned-char -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -std=c++11 -fno-strict-aliasing -fno-rtti

# disable flag
# -fno-exceptions
#for gcc 8+ add: -Wno-class-memaccess
GCCVER5 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 5)
ifeq "$(GCCVER5)" "1"
 BASEFLAGS += -Wno-implicit-fallthrough
endif

GCCVER8 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 8)
ifeq "$(GCCVER8)" "1"
  BASEFLAGS += -Wno-class-memaccess
endif

LINKER  := $(if $(LINKER),$(LINKER),g++)

LDFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)

# LDFLAGS += -L${BAM}

LIBS    := ${HTSLIB}/libhts.a ${LIBBZ2} ${LIBLZMA} ${LIBDEFLATE} -lz -lm

ifneq (,$(filter %nothreads %prof %profile, $(MAKECMDGOALS)))
 NOTHREADS=1
endif

#detect MinGW (Windows environment)
ifneq (,$(findstring mingw,$(shell ${CXX} -dumpmachine)))
 WINDOWS=1
endif

# Misc. system commands
#ifdef WINDOWS ##<-- use MSYS
# RM = del /Q
#else
RM = rm -f
#endif

# Non-windows systems need pthread
ifndef WINDOWS
 ifndef NOTHREADS
   LIBS := -pthread ${LIBS}
   BASEFLAGS += -pthread
 endif
endif

ifdef NOTHREADS
  BASEFLAGS += -DNOTHREADS
endif

# Compiling for Windows with MinGW?
#ifneq ($(findstring -mingw,$(shell $(CC) -dumpmachine 2>/dev/null)),)
#LIBS += -lregex -lws2_32
#endif
# File endings
ifdef WINDOWS
 EXE = .exe
 LIBS += -lregex -lws2_32
else
 EXE =
endif

DMACH := $(shell ${CXX} -dumpmachine)

ifneq (,$(filter %release %static %static-cpp, $(MAKECMDGOALS)))
  # -- release build
  RELEASE_BUILD=1
  CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O3)
  CXXFLAGS += -DNDEBUG $(BASEFLAGS)
else
  ifneq (,$(filter %memcheck %memdebug %tsan %tcheck %thrcheck, $(MAKECMDGOALS)))
     #use sanitizer in gcc 4.9+
     GCCVER49 := $(shell expr `${CXX} -dumpversion | cut -f1,2 -d.` \>= 4.9)
     ifeq "$(GCCVER49)" "0"
       $(error gcc version 4.9 or greater is required for this build target)
     endif
     CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
     SANLIBS :=
     ifneq (,$(filter %tsan %tcheck %thrcheck, $(MAKECMDGOALS)))
        # thread sanitizer only (incompatible with address sanitizer)
        CXXFLAGS += -fno-omit-frame-pointer -fsanitize=thread -fsanitize=undefined $(BASEFLAGS)
        SANLIBS := -ltsan
     else
        # address sanitizer
        CXXFLAGS += -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address $(BASEFLAGS)
        SANLIBS := -lasan
     endif
     ifeq "$(GCCVER5)" "1"
       CXXFLAGS += -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=vptr
       CXXFLAGS += -fsanitize=float-cast-overflow -fsanitize=object-size
       #CXXFLAGS += -fcheck-pointer-bounds -mmpx
     endif
     CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG -fno-common -fstack-protector
     LIBS := ${SANLIBS} -lubsan -ldl ${LIBS}
  else
     ifneq (,$(filter %prof %profile, $(MAKECMDGOALS)))
     ## profiling build
       CXXFLAGS := -DNDEBUG $(BASEFLAGS) -g -pg
       LDFLAGS += -g -pg
     else
        #just plain debug build
        DEBUG_BUILD=1
        CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
        ifneq (, $(findstring darwin, $(DMACH)))
           CXXFLAGS += -gdwarf-3
        endif
        CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
     endif
  endif
endif

ifdef RELEASE_BUILD
 ifneq ($(findstring static,$(MAKECMDGOALS)),) 
  # static or static-cpp found
  ifneq ($(findstring static-cpp,$(MAKECMDGOALS)),) 
     #not a full static build, only c/c++ libs
     LDFLAGS := -static-libgcc -static-libstdc++ ${LDFLAGS}
  else
     #full static build
     LDFLAGS := -static -static-libgcc -static-libstdc++ ${LDFLAGS}
  endif
 endif
endif

ifdef DEBUG_BUILD
  #$(warning Building DEBUG version of stringtie.. )
  DBG_WARN=@echo
  DBG_WARN+='WARNING: built DEBUG version [much slower], use "make clean release" for a faster, optimized version of the program.'
endif

OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o ${GDIR}/GSam.o \
 ${GDIR}/gdna.o ${GDIR}/codons.o ${GDIR}/GFastaIndex.o ${GDIR}/GFaSeqGet.o ${GDIR}/gff.o 


ifneq (,$(filter %memtrace %memusage %memuse, $(MAKECMDGOALS)))
    CXXFLAGS += -DGMEMTRACE
    OBJS += ${GDIR}/proc_mem.o
    OBJS_MULTI += ${GDIR}/proc_mem.o
endif

ifndef NOTHREADS
 OBJS += ${GDIR}/GThreads.o 
 OBJS_MULTI += ${GDIR}/GThreads.o 
endif

%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

# OBJS += rlink.o tablemaker.o t_record.o
OBJS += t_record.o dot_record.o rlink.o tablemaker.o unispg.o APPLY_UNISPG/unispg_A.o CREATE_UNISPG/unispg_C.o APPLY_UNISPG/infer_tranx_A.o CREATE_UNISPG/infer_tranx_C.o helper.o APPLY_UNISPG/parse_reads_A.o APPLY_UNISPG/processOptions_A.o CREATE_UNISPG/processOptions_C.o APPLY_UNISPG/processBundle_A.o CREATE_UNISPG/processBundle_C.o CREATE_UNISPG/build_graphs_C.o CREATE_UNISPG/rlink_C.o APPLY_UNISPG/rlink_A.o CREATE_UNISPG/rlink_C_help.o

all release static static-cpp debug: multistringtie${EXE}
memcheck memdebug tsan tcheck thrcheck: multistringtie${EXE}
memuse memusage memtrace: multistringtie${EXE}
prof profile: multistringtie${EXE}
nothreads: multistringtie${EXE}

multistringtie.o : $(GDIR)/GBitVec.h $(GDIR)/GHashMap.hh $(GDIR)/GSam.h

APPLY_UNISPG/processOptions_A.o: APPLY_UNISPG/processOptions_A.h definitions.h global_params.h
CREATE_UNISPG/processOptions_C.o: CREATE_UNISPG/processOptions_C.h definitions.h global_params.h

APPLY_UNISPG/processBundle_A.o: APPLY_UNISPG/processBundle_A.h definitions.h unispg.h rlink.h
CREATE_UNISPG/processBundle_C.o: CREATE_UNISPG/processBundle_C.h definitions.h unispg.h rlink.h


rlink.o : rlink.h tablemaker.h $(GDIR)/GSam.h $(GDIR)/GBitVec.h
CREATE_UNISPG/rlink_C.o : CREATE_UNISPG/rlink_C.h 
CREATE_UNISPG/rlink_C_help.o : CREATE_UNISPG/rlink_C_help.h 

APPLY_UNISPG/rlink_A.o : APPLY_UNISPG/rlink_A.h 

t_record.o : rlink.h t_record.h
tablemaker.o : tablemaker.h rlink.h
dot_record.o : tablemaker.h rlink.h

unispg.o : tablemaker.h
APPLY_UNISPG/unispg_A.o : tablemaker.h
CREATE_UNISPG/unispg_C.o : tablemaker.h


APPLY_UNISPG/infer_tranx_A.o: APPLY_UNISPG/infer_tranx_A.h rlink.h
CREATE_UNISPG/infer_tranx_C.o: CREATE_UNISPG/infer_tranx_C.h rlink.h

APPLY_UNISPG/parse_reads_A.o:
CREATE_UNISPG/build_graphs_C.o:

helper.o:

##${BAM}/libbam.a: 
##	cd ${BAM} && make lib

${HTSLIB}/libhts.a:
	cd ${HTSLIB} && ./build_lib.sh





multistringtie${EXE}: ${HTSLIB}/libhts.a $(OBJS) multistringtie.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}
test demo tests: multistringtie${EXE}
	@./run_tests.sh
.PHONY : clean cleanall cleanAll allclean

# target for removing all object files

#	echo $(PATH)
############################
## KH ADD
clean:
	${RM} multistringtie${EXE} multistringtie.o*  $(OBJS) $(OBJS_MULTI)

# pr: rm -rf multistringtie.o unispg.o
## END of KH ADD
############################
##allclean cleanAll cleanall:
##	cd ${BAM} && make clean
##	${RM} stringtie${EXE} stringtie.o* $(OBJS)
##	${RM} core.*
