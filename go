#gcc x.c -o x `pkg-config --cflags --libs gtk+-3.0`

#gcc -DRPF=1 x.c -o x `pkg-config --cflags --libs gtk+-2.0` -L/usr/lib
#gcc -DRPF=1 mm.c -o mm `pkg-config --cflags --libs gtk+-2.0` -L/usr/lib
#gcc -x c -DNATIVE=1 -DUNIX=1 -DGTK=1 -c -Wall -o alviewcore.o alviewcore.cpp

#gcc -c -Wall -o minibam.o minibam.c 

#gcc -c -Wall -o mit.o mit.c 
#gcc -x c -DUNIX=1 -DGTK=1 -c -Wall -I/data/nextgen/finneyr/samtools-0.1.18/ -o alviewcore.o alviewcore.cpp
#gcc -DRPF=1 alvgtk.c -o alvgtk.o `pkg-config --cflags --libs gtk+-2.0` -L/usr/lib

gcc -DUNIX=1 -DGTK=1 -g -x c -o alvgtk -I/home/rfinney/samtools-0.1.18/ alvgtk.c mit.c alviewcore.cpp `pkg-config --cflags --libs gtk+-2.0` -L/usr/lib -L/home/rfinney/samtools-0.1.18/ -lbam -lm -lz 



