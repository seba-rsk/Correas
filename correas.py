from math import pi, sqrt

class Material():
    def __init__(self):
        self.calidad = 'ZAR250'
        self.fy = 250            # Tensión de fluencia en MPa
        self.fu = 330            # Tensión de rotura en MPa
        self.mu = 0.30           # Módulo de Poisson en un período elástico [0.27-0.30]
        self.E = 200000          # Módulo de elásticidad longitudinal o módulo de Young en MPa
        self.G = 77200           # Módulo de elásticidad transversal en MPa
        self.gamma = 7.85        # Peso específico en t/m³

    def __str__(self):
        return f'Calidad del Acero: {self.calidad}\nTensión de fluencia fy = {self.fy} MPa\nTensión de rotura fu = {self.fu} MPa\nMódulo de Poisson µ = {self.mu}\nMódulo de elásticidad longitudinal E = {self.E} MPa\nMódulo de elásticidad transversal G = {self.G} MPa\nPeso específico γ = {self.gamma} t/m³\n'

    def change_material(self):
        pass

    def new_material(self,calidad,fy,fu):
        self.calidad = calidad
        self.fy = fy
        self.fu = fu

    def delete_material(self):
        pass


class Perfil_C(Material):
    def __init__(self, ht, bt, dt, t):
        Material.__init__(self)
        self.ht = ht    # Alto total en mm
        self.bt = bt    # Ancho total en mm
        self.dt = dt    # Labio total en mm
        self.t = t      # espesor en mm
        self.r = self.t    # Radio interior en mm
        self.hp = self.ht-2*self.t-2*self.r   # Alto plano alma en mm
        self.bp = self.bt-2*self.t-2*self.r   # Ancho plano ala en mm
        self.dtp = self.dt-self.t-self.r      # Alto plano ala en mm
        self.Ag = ((self.hp+self.bp*2+self.dtp*2)*self.t+pi*((self.r+self.t)**2-self.r**2))/100     # Sección transversal en cm²
        self.g = self.Ag*self.gamma/10        # Peso en kg/m
        self.Ix = (self.t*self.hp**3/12+2*((self.bp*self.t**3/12)+(self.t*self.bp)*(self.ht/2-self.t/2)**2)+2*((self.t*self.dtp**3/12)+(self.t*self.dtp)*(self.hp/2-self.dtp/2)**2)+2*pi/8*((self.t+self.r)**4-self.r**4)+pi*((self.t+self.r)**2-self.r**2)*(self.hp/2+(4/(3*pi))*((self.t+self.r)**3-self.r**3)/((self.t+self.r)**2-self.r**2))**2)/10000        # Módulo de Inercia respecto al eje X cm4
        self.Sx = self.Ix/(self.ht/20)      # Módulo Resistente respecto al eje X en cm³
        self.ix = sqrt(self.Ix/self.Ag)     # Radio de giro respecto al eje X en cm
        self.ey = (self.bt/2-((self.hp*self.t*(self.bt/2-self.t/2)-2*self.dtp*self.t*(self.bt/2-self.t/2))/(self.Ag*100)))/10       # Distancia al baricentro hozizontal en cm
        self.Iy = (self.hp*(self.t**3)/12+self.hp*self.t*(self.ey*10-self.t/2)**2+2*(self.t*self.bp**3/12+self.t*self.bp*(self.bt/2-self.ey*10)**2)+2*(self.dtp*self.t**3/12+self.dtp*self.t*(self.bt-self.t/2-self.ey*10)**2)+2*(pi/8-8/(9*pi))*((self.t+self.r)**4-self.r**4)+((self.ey*10-self.t-self.r+(4/(3*pi))*((self.r+self.t)**3-self.r**3)/(self.t**2+2*self.t*self.r))**2)*pi*(self.t**2+2*self.t*self.r)/2+((self.bt-self.ey*10-self.t-self.r+(4/(3*pi))*((self.r+self.t)**3-self.r**3)/(self.t**2+2*self.t*self.r))**2)*pi*(self.t**2+2*self.t*self.r)/2)/10000        # Módulo de Inercia respecto al eje Y en cm4
        self.Sy1 = self.Iy/self.ey                  # Módulo Resistente 1 respecto al eje Y en cm³
        self.Sy2 = self.Iy/(self.bt/10-self.ey)     # Módulo Resistente 2 respecto al eje Y en cm³
        self.iy = sqrt(self.Iy/self.Ag)             # Radio de giro respecto al eje Y en cm
        self.ec = (self.bt*(1+(2*self.dt/self.bt)*(1-(4*self.dt**2)/(3*self.ht**2)))/(2+self.ht/(3*self.bt)+(2*self.dt/self.bt)*(1-2*self.dt/self.ht+4*(self.dt**2)/(3*self.ht**2))))/10+self.ey-self.t/20       # Distancia del eje y al centro de corte en cm
        self.r0 = sqrt(self.ix**2+self.iy**2+self.ec**2)                # Radio de giro polar de la totalidad de la sección transversal respecto al centro de corte en cm
        self.Jt = ((2*self.dtp+2*self.bp+self.hp+pi*2*(self.r+0.5*self.t))*self.t**3)/30000       # Módulo de torsión en cm4
        self.Cw = ((((self.ht-self.t)**2)*((self.bt-self.t)**2)*self.t/12)*(2*((self.ht-self.t)**3)*(self.bt-self.t)+3*((self.ht-self.t)**2)*((self.bt-self.t)**2)+48*((self.dt-self.t/2)**4)+112*(self.bt-self.t)*((self.dt-self.t/2)**3)+8*(self.ht-self.t)*((self.dt-self.t/2)**3)+48*(self.ht-self.t)*(self.bt-self.t)*((self.dt-self.t/2)**2)+12*((self.ht-self.t)**2)*((self.dt-self.t/2)**2)+12*((self.ht-self.t)**2)*(self.bt-self.t)*(self.dt-self.t/2)+6*((self.ht-self.t)**3)*(self.dt-self.t/2))/(6*((self.ht-self.t)**2)*(self.bt-self.t)+((self.ht-self.t)+2*(self.dt-self.t/2))**3-24*(self.ht-self.t)*((self.dt-self.t/2)**2)))/1000000    # Módulo de alabeo en cm6


    def __str__(self):
        return f'Perfil C{self.ht}x{self.bt}x{self.dt}x{self.t}'

       
    def show_properties(self):
        print(f'''Propiedades geometricas:
              \t- ht = {round(self.ht,0)} mm
              \t- bt = {round(self.bt,0)} mm
              \t- dt = {round(self.dt,0)} mm
              \t- t = {round(self.t,2)} mm
              \t- r = {round(self.r,2)} mm
              \t- hp = {round(self.hp,2)} mm
              \t- bp = {round(self.bp,2)} mm
              \t- dtp = {round(self.dtp,2)} mm
              \t- Ag = {round(self.Ag,2)} cm²
              \t- g = {round(self.g,2)} kg/m
              \t- Ix = {round(self.Ix,2)} cm⁴
              \t- Sx = {round(self.Sx,2)} cm³
              \t- ix = {round(self.ix,2)} cm
              \t- ey = {round(self.ey,2)} cm
              \t- Iy = {round(self.Iy,2)} cm⁴
              \t- Sy1 = {round(self.Sy1,2)} cm³
              \t- Sy2 = {round(self.Sy2,2)} cm³
              \t- iy = {round(self.iy,2)} cm
              \t- ec = {round(self.ec,2)} cm
              \t- r0 = {round(self.r0,2)} cm
              \t- Jt = {round(self.Jt,2)} cm⁴
              \t- Cw = {round(self.Cw,2)} cm⁶''')


    def verif_geometrica(self):
     
        # Esbeltez labio (elemento no rigidizado - compresión con tensiones variables)
        self.esb_labio = self.dtp/self.t
        self.Is_labio = self.t*self.dtp**3/12 
        # Esbeltez Ala (elemento rigidizado por un ala y un labio - compresión con tensiones uniformes)
        self.esb_ala = self.bp/self.t
        self.Is_ala = self.t*self.bp**3/12 
        # Esbeltez Alma (elemento rigidizado - tensiones linealmente variables)
        self.esb_alma = self.hp/self.t
        if self.esb_labio < 60:         # B.1.1
            print(f'El labio verifica las consideraciones geometricas. dtp/t = {round(self.esb_labio,2)} < 60')
        else:
            print(f'El labio NO verifica las consideraciones geometricas. dtp/t = {round(self.esb_labio,2)} > 60')
        if self.esb_ala < 60:           # B.1.1
            print(f'El ala verifica las consideraciones geometricas. bp/t = {round(self.esb_ala,2)} < 60')
        else:
            print(f'El ala NO verifica las consideraciones geometricas. bp/t = {round(self.esb_ala,2)} > 60')
        if self.esb_alma < 200:         # B.1.2
            print(f'El alma verifica las consideraciones geometricas. hp/t = {round(self.esb_ala,2)} < 200')
        else:
            print(f'El labio NO verifica las consideraciones geometricas. hp/t = {round(self.esb_ala,2)} > 200')

    
    def seccion_efectiva(self):

        # Ancho efectivo labio (B.3.2 - Elemento no rigidizado con tensiones linealmente variables)
        self.Fcr_labio = 0.43*(pi**2*self.E)/(12*(1-self.mu**2)*(self.dtp/self.t)**2)       # Tensión de pandeo elástico de placas en el labio en MPa
        self.lamda_labio = sqrt(self.fy/self.Fcr_labio)
        self.ro_labio = (1-0.22/self.lamda_labio)/self.lamda_labio if self.lamda_labio>0.673 else 1
        self.dte = self.dtp*self.ro_labio                                             # Longitud efectiva del labio en mm
        
        # Ancho efectivo ala (B.4.2 - Elemento uniformemente comprimido con un rigidizador de borde)
        self.S = 1.28*sqrt(self.E/self.fy)
        if self.esb_ala <= 0.328*self.S:        # B.4.2-a
            self.be = self.bp                   # Longitud efectiuva del ala en mm
            self.be1 = self.be2 = self.be/2
            self.be_laguna = 0                  # Longitud no efectiva del ala en mm
            self.dte_reduc = self.dte
            self.dte_laguna = self.dtp-self.dte         # Longitud no efectiva del labio en mm
        else:
            self.Is_labio = self.t*self.dtp**3/12                                                                               # Inercia del rigidizador respecto del eje paralelo al elemento rigidizado en mm4
            self.Ia_labio = min(399*(self.t**4)*(self.esb_ala/self.S-0.328)**3, (self.t**4)*(115*self.esb_ala/self.S+5))        # Inercia necesario del rigidizador en mm4
            self.Ri = min(self.Is_labio/self.Ia_labio,1)
            self.n_ala = max(0.582-self.esb_ala/(4*self.S),1/3)
            self.D_b = self.dt/self.bp
            if self.D_b <= 0.25:
                self.k_ala = min(3.57*(self.Ri**self.n_ala)+0.43, 4)
            elif 0.25 < self.D_b <= 0.8:
                self.k_ala = min((4.82-5*self.D_b)*(self.Ri**self.n_ala)+0.43,4)
            else:
                self.k_ala = 4
                
            self.Fcr_ala = self.k_ala*(pi**2*self.E)/(12*(1-self.mu**2)*(self.bp/self.t)**2)        # Tensión de pandeo elástico de placas en el ala en MPa
            self.lamda_ala = sqrt(self.fy/self.Fcr_ala)
            self.ro_ala = (1-0.22/self.lamda_ala)/self.lamda_ala if self.lamda_ala>0.673 else 1
            self.be = self.bp*self.ro_ala               # Longitud efectiva del ala en mm
            self.be1 = (self.be/2)*self.Ri
            self.be2 = self.be-self.be1                                                                          
            self.be_laguna = self.bp-self.be            # Longitud no efectiva del ala en mm
            
            self.dte_reduc = self.dte*self.Ri
            self.dte_laguna = self.dtp-self.dte_reduc         # Longitud no efectiva del labio en mm

        # Ancho efectivo alma (B.2.3 - Elemento rigidizado con tensiones linealmente variables)
        self.psi = -1        # Ψ = abs(f1/f2) = 1 para flexión pura
        self.ht_bt = self.ht/self.bt
        self.beta_inicio = self.beta = 1            # Variables para iniciar la iteración
        for i in range(100000):                     # Bucle de iteración
            if self.beta > self.beta_inicio:        #Condición que controla la salida del bucle
                break
            else:
                self.beta_inicio = self.beta
           
            self.psi = -1+(i*1/100000)
            self.k_alma = 4+2*(1+abs(self.psi))**3+2*(1+abs(self.psi))                 # Coeficiente de abolladura
            self.Fcr_alma = self.k_alma*(pi**2*self.E)/(12*(1-self.mu**2)*(self.hp/self.t)**2)
            self.f1_alma = self.fy*(self.hp/2)/(self.ht/2)
            self.lamda_alma = sqrt(self.f1_alma/self.Fcr_alma)
            self.ro_alma = (1-0.22/self.lamda_alma)/self.lamda_alma if self.lamda_alma>0.673 else 1
            self.he = self.hp*self.ro_alma
            if self.ht_bt <= 4:
                self.he1 = self.he/(3+abs(self.psi))
                if abs(self.psi) > 0.236:
                    self.he2 = self.he/2
                else:
                    self.he2 = self.he-self.he1
            else:
                self.he1 = self.he/(3+abs(self.psi))
                self.he2 = self.he/(1+abs(self.psi))-self.he1
            self.he_laguna = max(self.hp/(1-min(self.psi,0))-self.he1-self.he2,0)
            # Sección efectiva
            self.Ae = self.Ag-self.t*(self.dte_laguna+self.be_laguna+self.he_laguna)/100        #Área efectiva en cm²
            self.yge = round(abs((self.t*(self.hp/2)*(self.hp/4)+self.bp*self.t*((self.ht/2)-self.t/2)+self.dtp*self.t*(self.ht/2-self.t-self.r-self.dtp/2)+pi*(self.t**2+2*self.t*self.r)*0.5*(self.ht/2-(self.r+self.t)+(4*((self.t+self.r)**3-self.r**3)/(3*pi*(self.t+self.r)**2-self.r*2)))-
                              self.t*(self.hp/2)*(self.hp/4)+self.t*self.he_laguna*(self.hp/2-self.he1-self.he_laguna/2)-self.bp*self.t*((self.ht/2)-self.t/2)+self.be_laguna*self.t*(self.ht/2-self.t/2)-self.dtp*self.t*(self.ht/2-self.t-self.r-self.dtp/2)+
                              self.dte_laguna*self.t*(self.hp/2-self.dtp+self.dte_laguna/2)-pi*(self.t**2+2*self.t*self.r)*0.5*(self.ht/2-(self.r+self.t)+(4*((self.t+self.r)**3-self.r**3)/(3*pi*(self.t+self.r)**2-self.r**2))))/(self.Ae*100)),4)       # [mm]
            self.yge1 = self.ht/2+self.yge
            self.yge2 = self.ht/2-self.yge
            self.alfa = -1*self.yge2/self.yge1
            self.beta = abs(self.alfa-self.psi)/abs(self.alfa)
        self.Jex = ((self.t*(self.hp/2+self.yge)**3)/12+self.t*(self.hp/2+self.yge)*((self.hp/2+self.yge)/2)**2-(self.t*self.he_laguna**3)/12-self.t*self.he_laguna*(self.hp/2+self.yge-self.he1-self.he_laguna/2)**2
                    +((self.bp-self.be_laguna)*self.t**3)/12+(self.bp-self.be_laguna)*self.t*(self.yge1-self.t/2)**2+self.t*((self.dtp-self.dte_laguna)**3)/12+self.t*(self.dtp-self.dte_laguna)*(self.yge1-self.t-self.r-(self.dtp-self.dte_laguna)/2)**2
                    +pi*0.125*((self.t+self.r)**4-self.r**4)+pi*((self.t+self.r)**2-self.r**2)*0.5*(self.yge1-(self.r+self.t)+(4*((self.t+self.r)**3-self.r**3)/(3*pi*((self.t+self.r)**2-self.r**2))))**2
                    +(self.t*(self.hp/2-self.yge)**3)/12+self.t*(self.hp/2-self.yge)*((self.hp/2-self.yge)/2)**2+(self.bp*self.t**3)/12+self.bp*self.t*(self.yge2-self.t/2)**2+self.t*(self.dtp**3)/12+self.t*self.dtp*(self.yge2-self.t-self.r-self.dtp/2)**2
                    +pi*0.125*((self.t+self.r)**4-self.r**4)+pi*((self.t+self.r)**2-self.r**2)*0.5*(self.yge2-(self.r+self.t)+(4*((self.t+self.r)**3-self.r**3)/(3*pi*((self.t+self.r)**2-self.r**2))))**2)/10000
        self.Sex=self.Jex/(self.yge1/10)

    def flexion_mayor_inercia(self, Mux, ly, lt, Cb):
        self.Mux = Mux
        self.ly = ly
        self.lt = lt
        self.Cb = Cb
        self.Lu = sqrt(0.36*self.Cb*pi**2*self.E*(self.ht/10)*self.Iy/(2*self.Sx*self.fy))                                                      # Long. de arriostramiento para la cual el PLT no es crítico.[cm]
        self.Sigma_ey = pi**2*self.E/(max(self.ly,0.0001)*100/self.iy)**2                                                                       # [MPa]
        self.Sigma_t = (self.G*self.Jt+pi**2*self.E*self.Cw/(max(self.lt,0.0001)*100)**2)/(self.Ag*self.r0**2)                                  # [MPa]
        self.Fe = self.Cb*self.r0*self.Ag*sqrt(self.Sigma_ey*self.Sigma_t)/self.Sx                                                              # [MPa]
        self.Fc = self.Fe if self.Fe<=0.56*self.fy else self.fy if self.Fe>=2.78*self.fy else (10/9)*self.fy*(1-10*self.fy/(36*self.Fe))        # [MPa]
        self.Mdx=0.95*self.fy*self.Sex/10000 if self.ly<self.Lu/100 else 0.9*self.Fc*self.Sex/10000
        self.ratio = self.Mux/self.Mdx

    def show_seccion_efectiva(self):
        print('Ancho efectivo labio')
        print(f'Fcr labio = {round(self.Fcr_labio,2)} MPa')
        print(f'λ labio = {round(self.lamda_labio,2)}')
        print(f'ρ labio = {round(self.ro_labio,2)}')
        print(f'dt efectivo = {round(self.dte,2)} mm')
        print(f'\nAncho efectivo ala')
        print(f'S = {round(self.S,2)}')
        if self.esb_ala <= 0.328*self.S:        # B.4.2-a
            print(f'be = {round(self.be,2)} mm')
            print(f'be1 = be2 = {round(self.be2,2)} mm')
            print(f'be_laguna = {round(self.be_laguna,2)} mm')
            print(f'dte_reduc = {round(self.dte_reduc,2)} mm')
            print(f'dte_laguna = {round(self.dte_laguna,2)} mm')
        else:
            print(f'Is labio = {round(self.Is_labio,2)} mm⁴')
            print(f'Ia labio = {round(self.Ia_labio,2)} mm⁴')
            print(f'Ri = {round(self.Ri,2)}')
            print(f'n ala = {round(self.n_ala,2)}')
            print(f'D/b = {round(self.D_b,2)}')
            if self.D_b <= 0.25:
                print(f'k ala = {round(self.k_ala,2)}')   
            elif 0.25 < self.D_b <= 0.8:
                print(f'k ala = {round(self.k_ala,2)}')
            else:
                print(f'k ala = {round(self.k_ala,2)}')
            print(f'Fcr ala = {round(self.Fcr_ala,2)} MPa')
            print(f'λ ala = {round(self.lamda_ala,2)}')
            print(f'ρ ala = {round(self.ro_ala,2)}')
            print(f'be = {round(self.be,2)} mm')
            print(f'be1 = {round(self.be1,2)} mm')
            print(f'be2 = {round(self.be2,2)} mm')                                                                     
            print(f'be_laguna = {round(self.be_laguna,2)} mm')
            print(f'dte reducido = {round(self.dte_reduc,2)} mm')
            print(f'dte_laguna = {round(self.dte_laguna,2)} mm')
        print(f'\nAncho efectivo alma')        
        print(f'Ψ = {round(self.psi,2)}')
        print(f'k alma = {round(self.k_alma,2)}')
        print(f'Fcr alma = {round(self.Fcr_alma,2)} MPa')
        print(f'f1 alma = {round(self.f1_alma,2)} Mpa')
        print(f'λ alma = {round(self.lamda_alma,2)}')
        print(f'ρ alma = {round(self.ro_alma,2)}')
        print(f'he = {round(self.he,2)} mm')
        print(f'ht/bt = {round(self.ht_bt,2)}')
        print(f'he1 = {round(self.he1,2)} mm') 
        print(f'he2 = {round(self.he2,2)} mm')
        print(f'he_laguna = {round(self.he_laguna,2)} mm')
        print(f'\nSección efectiva') 
        print(f'Ae = {round(self.Ae,2)} cm²')
        print(f'yge = {round(self.yge,2)} mm')
        print(f'yge1 = {round(self.yge1,2)} mm')
        print(f'yge2 = {round(self.yge2,2)} mm')
        print(f'Jex = {round(self.Jex,2)} cm⁴')
        print(f'Sex = {round(self.Sex,2)} cm³')


    def show_flexion_mayor_inercia(self):
        print('Input:')
        print(f'Mux = {round(self.Mux,2)} tm')
        print(f'ky.Ly = {round(self.ly,2)} m')
        print(f'kt.Lt = {round(self.lt,2)} m')
        print(f'Cb = {round(self.Cb,2)}')
        print(f'\nVerificación:')
        print(f'Lu = {round(self.Lu,2)} cm')
        print(f'Sigma_ey = {round(self.Sigma_ey,2)} MPa')
        print(f'Sigma_t = {round(self.Sigma_t,2)} MPa')
        print(f'Fe = {round(self.Fe,2)} MPa')
        print(f'Fc = {round(self.Fc,2)} MPa')
        print(f'Mdx = {round(self.Mdx,2)} Mpa')
        print(f'ratio = {round(self.ratio,2)}')
