import math
import pylab

print('如果不输入值，则自动设为默认值')
todefault = input('是否直接使用所有默认值？（Y/N)')


def inputvalue(word, default='', de=todefault):
    if de == 'Y':
        return default
    else:
        input_value = input(word)
        if input_value == '':
            return default
        else:
            return input_value


def roughness(f, aoi):
    AoG = 3.1416 / 2 - aoi
    a = (0.006 + 0.0102 * VoW) ** (0-1)
    v31 = math.sin(AoG) * (1 - (math.exp(-a * (AoG ** 2) / 4)) / (AoG * ((3.1416 * a) ** 0.5)))
    v32 = math.sin(AoG) / 2
    v3 = max(v31, v32)
    SL = -10 * math.log10(1 - v3) - 20 * math.log10(0.3 + 0.7 / (1 + 6 * (10 ** -11) * (VoW ** 4) * (f ** 2)))
    CoR = 10 ** (-SL / 20)
    return CoR


def pluslist(a,b):
    try:
        t = []
        for i in range(len(a)):
            t.append(a[i]+b[i])
    except:
        t = a
    return t


Rho = 1024
c = 1500
NoF = 100
NoG = 25
P = 1
try:
    ds = float(inputvalue('声源深度（单位：米，默认4）：', '4'))
    H = float(inputvalue('水深（单位：米，默认60）：', '60'))
    dCPA = float(inputvalue('测量距离（单位：米，默认200）：', '200'))
    NoH = int(inputvalue('水听器数量（默认3）：', '3'))
    HTop = float(inputvalue('顶部水听器深度（单位：米，默认15）：', '15'))
    DeltaH = float(inputvalue('水听器间距（单位：米，默认20）：', '20'))
    VoW = float(inputvalue('风速（单位：节，默认0）：', '0'))
    material = inputvalue('海底材质（默认细沙）：', '细沙')
    if material == '细沙':
        RhoSB = 1941
        cSB = 1749
        Att = 0.9
    else:
        print('未收录该材料,请输入参数：')
        RhoSB = float(input('密度(单位kg/m^3)：'))
        cSB = float(input('声速(仅支持>1500m/s的材料）：'))
        Att = float(input('衰减(单位dB/波长）：'))
except:
    print('非法输入，请重新运行')
    exit()
CGA = 3.1416/2 - math.asin(c / cSB)
f13 = [10**(19/20+x/10) for x in range(51)]
f = [10**(x/10/NoF+(1/NoF + 19)/20) for x in range(50*NoF)]
k = [2*3.1416*x/c for x in f]
dH = [-HTop-DeltaH*x for x in range(NoH)]
r = [((ds+dH[x])**2+dCPA**2)**0.5 for x in range(NoH)]


def reflection(aoi, c_water=c, r_water=Rho, c_seabed=cSB, r_seabed=RhoSB, attenuation=Att):
    alphaC = attenuation / 8.686 / c_seabed
    c_complex = 2 * 3.1416 / (2 * 3.1416 / cSB - 1j * alphaC)
    sinaoic = c_complex * math.sin(aoi) / c_water
    cosaoic = (1-sinaoic**2)**0.5
    z1 = c_water * r_water
    z2 = c_complex * r_seabed
    coefficient = (z2 * math.cos(aoi) - z1 * cosaoic) / (z2 * math.cos(aoi) + z1 * cosaoic)
    return coefficient


def pressure(noh, nog, fcount):
    pis = [-2*nog*H-ds, 2*nog*H+ds, 2*(nog+1)*H-ds, -2*(nog+1)*H+ds]
    rlocal = [((pis[x]-dH[noh])**2 + dCPA**2)**0.5 for x in range(4)]
    aoilocal = [math.asin(dCPA/x) for x in rlocal]
    prz = P*(((-reflection(aoilocal[0]) * roughness(f[fcount], aoilocal[0])) ** nog)
           * (2.71828**(0 - k[fcount] * rlocal[0] * 1j)) / rlocal[0] -
           ((-reflection(aoilocal[1]) * roughness(f[fcount], aoilocal[1])) ** nog) * roughness(f[fcount], aoilocal[1])
           * (2.71828**(0 - k[fcount] * rlocal[1] * 1j)) / rlocal[1] +
           ((-reflection(aoilocal[2]) * roughness(f[fcount], aoilocal[2])) ** (nog+1))
           * (2.71828**(0 - k[fcount] * rlocal[2] * 1j)) / rlocal[2] -
           ((-reflection(aoilocal[3]) * roughness(f[fcount], aoilocal[3])) ** (nog+1)) / roughness(f[fcount], aoilocal[3])
           * (2.71828**(0 - k[fcount] * rlocal[3] * 1j)) / rlocal[3])
    return prz


p4h = []
for i in range(NoH):
    sump = []
    for ii in range(len(f)):
        sump.append(((r[i] * abs(sum([pressure(i, x, ii) for x in range(NoG)])) / P)**2))
    p4h.append(sump)
ios = []
for i in range(NoH):
    poeh = []
    for ii in range(len(f13)-1):
        iii = ii * NoF
        poeh.append((sum(p4h[i][iii:iii + NoF])/NoF))
    ios = pluslist(poeh, ios)
x=[]
y=[]
for i in range(len(ios)):
    y.append(10*math.log10(ios[i]/NoH))
    x.append((f13[i]*f13[i+1])**0.5)
pylab.plot(x,y)
pylab.xlabel(u"frequency")
pylab.ylabel(u"correction factor")
pylab.title(u'delta-f curve')
pylab.semilogx()
file = open(r'./record.txt', 'w')
record=['x, y']
for i in range(len(x)):
    record.append(str(x[i])+' '+str(y[i]))
file.writelines(record)
pylab.show()
