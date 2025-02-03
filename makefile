ADCPREP = /workspace/adcirc/adc/adcprep 
PADCIRC = /workspace/adcirc/adc/padcirc
ADCPREP_DG = /workspace/adcirc/build/adcprep
PDG = /workspace/adcirc/build/padcirc
DEBUG = /workspace/adcirc/debug/adcirc

NP = 8

shinnecock: shinnecock-dg shinnecock-adc

shinnecock-debug:
	cd shinnecock && \
	gdb $(DEBUG)

shinnecock-serial:
	cd shinnecock && \
	$(DEBUG)

shinnecock-dg:
	cd shinnecock && \
	rm -rf PE* && \
	$(ADCPREP_DG) --np $(NP) --partmesh && \
	$(ADCPREP_DG) --np $(NP) --prepall && \
	mpirun --np $(NP) $(PDG) && \
	mv fort.61 fort.61.adg

shinnecock-adc:
	cd shinnecock && \
	rm -rf PE* && \
	$(ADCPREP) --np $(NP) --partmesh && \
	$(ADCPREP) --np $(NP) --prepall && \
	mpirun --np $(NP) $(PADCIRC) && \
	mv fort.61 fort.61.adc
