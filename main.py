import correas

# Programa principal
material_1 = correas.Material()
print('MATERIAL SELECCIONADO:')
print(material_1)
print('\n')

# H [mm], B [mm], D [mm], t [mm]
perfil = correas.Perfil_C(400, 120, 120, 2.0)
print('PERFIL SELECCIONADO:')
print(perfil)
print('\n')
print('PROPIEDADES GEOMETRICAS:')
perfil.show_properties()
print('\n')
print('VERIFICACION GEOMETRICA')
perfil.verif_geometrica()
print('\n')
print('SECCION EFECTIVA')
perfil.seccion_efectiva()
perfil.show_seccion_efectiva()
print('\n')
print('VERIFICACION A FLEXION')

# Mux [tm], ky.Ly [m], kt.Lt [m], Cb
perfil.flexion_mayor_inercia(2, 3, 3, 1)
perfil.show_flexion_mayor_inercia()
