Definition of the event topologies
-----------------------------------------------
-----------------------------------------------

The topologies for the categorization of the events are defined in ```/src/SignalDefinition.cc```. 
Starting from topology = 0, a bit is set for every different category.

The first bit (0) is set for all the events with the right ZZ &#8594 4l Z bosons (correct mass range and right opposite sign-same flavour leptons pairing). For all the events out of this category (with the first bit = 0) the other bits are not set.

In the following the position of the bit set (left) and the event category (right) are listed:

- 0 &#8594 ZZ4l 

- 1 &#8594 + jets(Right Pt and eta range)
  
- 2 &#8594 + 2q
  
- 3 &#8594 + 1lepton
  
- 4 &#8594 + hadronic W
    
- 5 &#8594 + hadronic Z
