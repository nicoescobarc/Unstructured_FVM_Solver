# Asigna el nombre del archivo de salida a la variable TARGET
TARGET = ConvDifReac

# Asigna el comando que llama al compilador (en este caso gfortran) a la variable "CC"
CC = gfortran

# Define las 'opciones' que se le pasan al compilador (flags)
CFLAGS = -O3

# Define los 'objetos' (o archivos) necesarios para compilar el objetivo (o 'target') principal
objects = constants_module.o subroutines_module.o $(TARGET).o


# Este es el objetivo principal del makefile. Por defecto es el unico objetivo que se ejecuta.
# Los demas objetivos sólo se ejecutan si estos son prerrequisitos del objetivo principal, en
# este caso esto se da cuando pertenecen a la variable "objects".
$(TARGET).exe: $(objects)
	$(CC) $(CFLAGS) -o $(TARGET).exe $(objects)

# Esta es una regla resume la aplicación de la misma regla para todos los objetos de extensión ".o" (que son los prerrequisitos del objetivo principal); es decir, para los que estén en la variable "objects".
# La regla dice que cada archivo ".o" requerido por el objetivo principal depende de un archivo "f.95" que tiene el mismo nombre, luego dice que para generar el archivo ".o" se compila el archivo ".f95" con el compilador definido en "CC" y los flags definidos en "CFLAGS", el flag adicional "-o" indica que el nombre del archivo ".o" resultante se especifica en la variable "$@"; esta variable le indica que el archivo resultante se nombrará de acuerdo al nombre de la variable a la izquierda de los dos puntos (del ":").
# Finalmente el flag "-c" indica que no se generará un ejecutable sino que se va a generar el archivo objeto (".o"), usando el archivo ".f95" que le indica la variable "$<"; esta variable toma el nombre del primer item en la lista de prerrequisitos, nótese que en este caso solo hay un prerrequisito (que es el archivo ".f95") para generar el archivo ".o" (es decir: para generar a YYY.o solo se necesita a YYY.f95)
%.o: %.f95
	$(CC) $(CFLAGS) -o $@ -c $<

# Esto indica que al teclear make clean se ejecutará ese comando en la consola
clean:
	rm -rf *.mod *.o *.exe










