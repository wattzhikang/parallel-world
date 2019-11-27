all : 
	cd src && $(MAKE) parallel-world

.PHONY : doc

doc : 
	cd src && $(MAKE) doc
