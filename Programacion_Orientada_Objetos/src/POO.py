'''
NAME
    POO

VERSION
    1.0        

AUTHOR
    García González Carlos        

DESCRIPTION
    Este código muestra cómo crear una jerarquía de clases con herencia, 
    cómo sobrescribir métodos y cómo utilizar un diccionario para almacenar y 
    acceder a objetos.

CATEGORY
    --        

USAGE

    % python programName [-parameter] [parameter_value]
    

ARGUMENTS
    No tiene metodos

METHOD


SEE ALSO


        
'''


# ===========================================================================
# =                            imports
# ===========================================================================





# ===========================================================================
# =                            Command Line Options
# ===========================================================================






# ===========================================================================
# =                            functions
# ===========================================================================




# ===========================================================================
# =                            main
# ===========================================================================

# Definimos la clase base Animal
class Animal:
    # Constructor de la clase Animal con atributos nombre y edad
    def __init__(self, nombre, edad):
        self.nombre = nombre
        self.edad = edad

    # Método haz_ruido que imprime un ruido genérico
    def haz_ruido(self):
        print("AAAAAAAAAAAAH")

# Definimos la clase Perro que hereda de Animal
class Perro(Animal):
    # Constructor de la clase Perro que añade el atributo raza y llama al constructor de la clase base
    def __init__(self, nombre, edad, raza):
        super().__init__(nombre, edad)
        self.raza = raza

    # Sobrescribimos el método haz_ruido para que imprima "Guau"
    def haz_ruido(self):
        print("Guau")

# Definimos la clase Gato que hereda de Animal
class Gato(Animal):
    # Constructor de la clase Gato que añade el atributo usa_arenero y llama al constructor de la clase base
    def __init__(self, nombre, edad, usa_arenero):
        super().__init__(nombre, edad)
        self.usa_arenero = usa_arenero

    # Sobrescribimos el método haz_ruido para que imprima "Miau"
    def haz_ruido(self):
        print("Miau")

# Crear objetos de las clases Perro y Gato
mi_perro = Perro("Gidget", 3, "Poddle")
mi_gato = Gato("Clementina", 2, True)

# Llamar a las funciones de los objetos y mostrar sus atributos
print(f"Perro: {mi_perro.nombre}, Edad: {mi_perro.edad}, Raza: {mi_perro.raza}")
mi_perro.haz_ruido()  # Esto imprimirá "Guau"

print(f"Gato: {mi_gato.nombre}, Edad: {mi_gato.edad}, Usa arenero: {'Sí' if mi_gato.usa_arenero else 'No'}")
mi_gato.haz_ruido()  # Esto imprimirá "Miau"

# Usar un diccionario para almacenar y acceder a los objetos
animales_dict = {
    'perro': mi_perro,
    'gato': mi_gato
}

# Acceder a los objetos en el diccionario y llamar a sus métodos
print("\nAccediendo a través del diccionario:")
for tipo, animal in animales_dict.items():
    print(f"Tipo: {tipo}, Nombre: {animal.nombre}, Edad: {animal.edad}")
    animal.haz_ruido()
