all: teste

teste: main.c funcs.c funcs.h
	cc -I/opt/local/include main.c funcs.c -lm -L/opt/local/lib -lnetpbm -o teste

clean:
	rm -f teste
