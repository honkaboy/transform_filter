PROJECT = transform_filter

SOURCES = \
  main.c

OBJECTS = $(SOURCES:.c=.o)

.phony: all
all: $(PROJECT).exe 

.phony: clean
clean:
	@rm -f $(PROJECT).bin $(OBJECTS)

$(PROJECT).exe: $(OBJECTS)
	@g++ $(OBJECTS) -o $@

