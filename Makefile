# -*- Makefile -*-
SHELL=/bin/sh

default:
	$(MAKE) $(MFLAGS) -C D1
	$(MAKE) $(MFLAGS) -C D2
	# $(MAKE) $(MFLAGS) -C D3

check:
	$(MAKE) $(MFLAGS) -C D1 check
	$(MAKE) $(MFLAGS) -C D2 check
	# $(MAKE) $(MFLAGS) -C D3 check

clean:
	$(MAKE) $(MFLAGS) -C D1 clean
	$(MAKE) $(MFLAGS) -C D2 clean
	# $(MAKE) $(MFLAGS) -C D3 clean
