PROGRAM = irvaudit

RM = rm -rf
OBJDIR = obj

BASEDIRS = -I. 


INCLUDEDIRS = $(BASEDIRS)

CXX = g++
LD =
SUFFIX = o

CXXFLAGS = -Wall -std=c++11 -pedantic -g $(INCLUDEDIRS) -m64 -fPIC \
	-fexceptions -DNEBUG -DIL_STD -Wno-long-long \
	-Wno-attributes  -fpermissive -Wno-sign-compare


LDFLAGS =  -lboost_system  -lboost_filesystem -lrt 

RENAME = -o

CXXSOURCES = \
	irvaudit.cpp \
	sim_irv.cpp \
	model.cpp  \
	audit.cpp 
	
CXXOBJECTS = $(patsubst %.cpp, $(OBJDIR)/%.$(SUFFIX), $(CXXSOURCES))

all : $(PROGRAM)

$(PROGRAM) : $(CXXOBJECTS)
	$(CXX) -o ${@} $(CXXOBJECTS) $(LD) $(LDFLAGS) 

$(OBJDIR)/%.$(SUFFIX) : %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(RENAME) $(@D)/$(@F) -c $(<)

clean:
	$(RM) $(CXXOBJECTS) $(PROGRAM) $(OBJDIR)


