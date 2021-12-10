all: teste

teste: main.c funcs.c funcs.h 
	cc -I/opt/local/include main.c funcs.c -lm -L/opt/local/lib -lnetpbm -o teste

teste2: main.c funcs.c funcs.h 
	cc -I/opt/local/include main.c funcs.c -L/usr/include/netpbm -o teste2


clean:
	rm -f teste
