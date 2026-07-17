include Makefile

SVN_TAG = $(shell svn info |grep 'Last Changed Rev' | awk -F: '{print $$2}' | sed -e 's/ *//g')
PACKDIR = 4server-both-$(SVN_TAG)

pack4server:
	svn export . ../$(PACKDIR)
	mkdir ../$(PACKDIR)/output
	cp fesigma.nc fe.nc foj.nc ftj.nc ../$(PACKDIR)
	cd ..; tar -czf $(PACKDIR).tgz $(PACKDIR)/*

test4server: 
	cd ..; tar -zxf $(PACKDIR).tgz
	cd ../$(PACKDIR); make DBG=1 frun
	cd ../$(PACKDIR); make DBG=1 run

remove:
	rm -r ../$(PACKDIR)
	rm -r ../$(PACKDIR).tgz
