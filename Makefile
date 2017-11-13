DIR=src
LIBS=-lm
OBJS_POP= $(patsubst %.c,%.o,$(DIR)/interface.c) $(patsubst %.c,%.o,$(DIR)/mylib.c) $(patsubst %.c,%.o,$(DIR)/pop_abc.c)\
 $(patsubst %.c,%.o,$(DIR)/pop_convertabc2.c) $(patsubst %.c,%.o,$(DIR)/pop_convertabc3.c) $(patsubst %.c,%.o,$(DIR)/pop_convertabc4.c)\
 $(patsubst %.c,%.o,$(DIR)/pop_convertabc.c) $(patsubst %.c,%.o,$(DIR)/pop_firstpass.c) $(patsubst %.c,%.o,$(DIR)/pop_genetictree.c)\
 $(patsubst %.c,%.o,$(DIR)/pop_joindata.c) $(patsubst %.c,%.o,$(DIR)/pop_makepop.c) $(patsubst %.c,%.o,$(DIR)/pop_makeprior.c)\
 $(patsubst %.c,%.o,$(DIR)/pop_makestats.c) $(patsubst %.c,%.o,$(DIR)/pop_maketarget.c) $(patsubst %.c,%.o,$(DIR)/pop_samplepriors.c)\
 $(patsubst %.c,%.o,$(DIR)/pop_summstats.c) 
OBJS_REJ=$(patsubst %.c,%.o,$(DIR)/firstpass.c) $(patsubst %.c,%.o,$(DIR)/mylib.c)
OBJS_SHU=$(patsubst %.c,%.o,$(DIR)/mylib.c) $(patsubst %.c,%.o,$(DIR)/turnit.c)
OBJS_SIM=$(patsubst %.c,%.o,$(DIR)/abc.c) $(patsubst %.c,%.o,$(DIR)/genetictree.c) $(patsubst %.c,%.o,$(DIR)/mylib.c)\
 $(patsubst %.c,%.o,$(DIR)/samplepriors.c) $(patsubst %.c,%.o,$(DIR)/summstats.c)
OBJS_SUM=$(patsubst %.c,%.o,$(DIR)/maketarget.c) $(patsubst %.c,%.o,$(DIR)/mylib.c) $(patsubst %.c,%.o,$(DIR)/summstats.c)

all: popabc rejection shuffle simulate summdata

popabc: $(OBJS_POP)
	$(CXX) $(OBJS_POP) -o $@.exe $(LIBS)

rejection: $(OBJS_REJ)
	$(CXX) $(OBJS_REJ) -o $@.exe $(LIBS)

shuffle: $(OBJS_SHU)
	$(CXX) $(OBJS_SHU) -o $@.exe $(LIBS)

simulate: $(OBJS_SIM)
	$(CXX) $(OBJS_SIM) -o $@.exe $(LIBS)

summdata: $(OBJS_SUM)
	$(CXX) $(OBJS_SUM) -o $@.exe $(LIBS)

clean:
	rm $(wildcard $(DIR)/*.o) popabc.exe rejection.exe shuffle.exe simulate.exe summdata.exe
