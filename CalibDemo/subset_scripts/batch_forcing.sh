#!/bin/bash
for i in 2007 2008 2009 2010 2011 2012 2013 2014 2015; do
   cp script_forcing_subset.txt script_forcing_subset${i}.txt
   sed -i "2s/.*/OLDFORCPATH=\'\/glade\/scratch\/zhangyx\/WRF-Hydro\/ForcingData.symlinks\/2007-2016\/${i}*LDASIN*\'/" script_forcing_subset${i}.txt
   nohup ./script_forcing_subset${i}.txt > nohup${i}.out 2>&1&
done
