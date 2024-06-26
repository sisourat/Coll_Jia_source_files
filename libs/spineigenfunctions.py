import sys
import numpy as np

def printCSF(spin,paired,unpaired,irun):

# N=0 CLOSED SHELL SINGLET
    if(spin==1 and len(unpaired)==0):
        ireturn = ne0spin1(paired,unpaired,irun)
# N=1 DOUBLET
    elif(spin==2 and len(unpaired)==1):
        ireturn = ne1spin2(paired,unpaired,irun)
# N=2 SINGLET
    elif(spin==1 and len(unpaired)==2):
        ireturn = ne2spin1(paired,unpaired,irun)
# N=2 TRIPLET
    elif(spin==3 and len(unpaired)==2):
        ireturn = ne2spin3(paired,unpaired,irun)
# N=3 QUADRUPLET
    elif(spin==4 and len(unpaired)==3):
        ireturn = ne3spin4(paired,unpaired,irun)
# N=3 DOUBLET
    elif(spin==2 and len(unpaired)==3):
        ireturn = ne3spin2(paired,unpaired,irun)
# N=4 SINGLET
    elif(spin==1 and len(unpaired)==4):
        ireturn = ne4spin1(paired,unpaired,irun)
# N=4 TRIPLET
    elif(spin==3 and len(unpaired)==4):
        ireturn = ne4spin3(paired,unpaired,irun)
# N=5 DOUBLET
    elif(spin==2 and len(unpaired)==5):
        ireturn = ne5spin2(paired,unpaired,irun)
    else:
        ireturn = 0
#        print("CSFs Not Implemented, I stop")
#        print(spin, len(unpaired)
#        sys.exit()

    return ireturn

def ne0spin1(paired,unpaired,irun):
    flist = open('list'+str(irun)+'.txt','a')
    npaired = len(paired)
    nunpaired = len(unpaired)
    if(not(nunpaired==0)):
        print("error in ne0spin1, I stop")
        sys.exit()

    print(1, file=flist)
    print(1.0," a ",'  '.join(map(str, paired)),'  '.join(map(str, unpaired))," b ",'  '.join(map(str, paired)), file=flist)
    return 1

def ne1spin2(paired,unpaired,irun):
    flist = open('list'+str(irun)+'.txt','a')
    npaired = len(paired)
    nunpaired = len(unpaired)
    if(not(nunpaired==1)):
        print("error in ne1spin2, I stop")
        sys.exit()

    print(1, file=flist)
    print(1.0," a ",'  '.join(map(str, paired)),'  '.join(map(str, unpaired))," b ",'  '.join(map(str, paired)), file=flist)
    return 1

def ne2spin1(paired,unpaired,irun):
    flist = open('list'+str(irun)+'.txt','a')
    npaired = len(paired)
    nunpaired = len(unpaired)
    if(not(nunpaired==2)):
        print("error in ne2spin1, I stop")
        sys.exit()

    print(2, file=flist)
    print("0.70710678118654746"," a ",'  '.join(map(str, paired)),str(unpaired[0])," b ",'  '.join(map(str, paired)), str(unpaired[1]), file=flist)
    print("0.70710678118654746"," a ",'  '.join(map(str, paired)),str(unpaired[1])," b ",'  '.join(map(str, paired)), str(unpaired[0]), file=flist)
    return 1

def ne2spin3(paired,unpaired,irun):
    flist = open('list'+str(irun)+'.txt','a')
    npaired = len(paired)
    nunpaired = len(unpaired)
    if(not(nunpaired==2)):
        print("error in ne2spin3, I stop")
        sys.exit()

    print(1, file=flist)
    print(1.0," a ",'  '.join(map(str, paired)),'  '.join(map(str, unpaired))," b ",'  '.join(map(str, paired)), file=flist)
    return 1

def ne3spin4(paired,unpaired,irun):
    flist = open('list'+str(irun)+'.txt','a')
    npaired = len(paired)
    nunpaired = len(unpaired)
    if(not(nunpaired==3)):
        print("error in ne3spin4, I stop")
        sys.exit()

    print(1, file=flist)
    print(1.0," a ",'  '.join(map(str, paired)),'  '.join(map(str, unpaired))," b ",'  '.join(map(str, paired)), file=flist)
    return 1

def ne3spin2(paired,unpaired,irun):
    flist = open('list'+str(irun)+'.txt','a')
    npaired = len(paired)
    nunpaired = len(unpaired)
    if(not(nunpaired==3)):
        print("error in ne3spin2, I stop")
        sys.exit()

    print(2, file=flist)
    print("0.70710678118654746"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[2])," b ",'  '.join(map(str, paired)),str(unpaired[1]), file=flist)
    print("0.70710678118654746"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[2])," b ",'  '.join(map(str, paired)),str(unpaired[0]), file=flist)
    print(3, file=flist)
    print("-0.40824829046386307"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[2])," b ",'  '.join(map(str, paired)),str(unpaired[0]), file=flist)
    print("0.40824829046386307"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[2])," b ",'  '.join(map(str, paired)),str(unpaired[1]), file=flist)
    print("0.81649658092772615"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[1])," b ",'  '.join(map(str, paired)),str(unpaired[2]), file=flist)
    return 2

def ne4spin1(paired,unpaired,irun):
    flist = open('list'+str(irun)+'.txt','a')
    npaired = len(paired)
    nunpaired = len(unpaired)
    if(not(nunpaired==4)):
        print("error in ne4spin1, I stop")
        sys.exit()

    print(4, file=flist)
    print(0.5," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[2])," b ",'  '.join(map(str, paired)), str(unpaired[1]), str(unpaired[3]), file=flist)
    print(0.5," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[3])," b ",'  '.join(map(str, paired)), str(unpaired[1]), str(unpaired[2]), file=flist)
    print(0.5," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[2])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[3]), file=flist)
    print(0.5," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[3])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[2]), file=flist)
    print(6, file=flist)
    print("0.57735026918962584"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[1])," b ",'  '.join(map(str, paired)), str(unpaired[2]), str(unpaired[3]), file=flist)
    print("0.57735026918962584"," a ",'  '.join(map(str, paired)),str(unpaired[2]),str(unpaired[3])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[1]), file=flist)
    print("0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[2])," b ",'  '.join(map(str, paired)), str(unpaired[1]), str(unpaired[3]), file=flist)
    print("0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[3])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[2]), file=flist)
    print("-0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[3])," b ",'  '.join(map(str, paired)), str(unpaired[1]), str(unpaired[2]), file=flist)
    print("-0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[2])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[3]), file=flist)
    return 2

def ne4spin3(paired,unpaired,irun):
    flist = open('list'+str(irun)+'.txt','a')
    npaired = len(paired)
    nunpaired = len(unpaired)
    if(not(nunpaired==4)):
        print("error in ne4spin1, I stop")
        sys.exit()

    print(2, file=flist)
    print("0.70710678118654746"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[2]),str(unpaired[3])," b ",'  '.join(map(str, paired)),str(unpaired[1]), file=flist)
    print("0.70710678118654746"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[2]),str(unpaired[3])," b ",'  '.join(map(str, paired)),str(unpaired[0])	, file=flist)
    print(3, file=flist)
    print("-0.40824829046386307"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[2]),str(unpaired[3])," b ",'  '.join(map(str, paired)),str(unpaired[0]), file=flist)
    print("0.40824829046386307"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[2]),str(unpaired[3])," b ",'  '.join(map(str, paired)),str(unpaired[1]), file=flist)
    print("0.81649658092772615"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[1]),str(unpaired[3])," b ",'  '.join(map(str, paired)),str(unpaired[2]), file=flist)
    print(4, file=flist)
    print("-0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[2]),str(unpaired[3])," b ",'  '.join(map(str, paired)),str(unpaired[0]), file=flist)
    print("0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[2]),str(unpaired[3])," b ",'  '.join(map(str, paired)),str(unpaired[1]), file=flist)
    print("-0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[1]),str(unpaired[3])," b ",'  '.join(map(str, paired)),str(unpaired[2]), file=flist)
    print("-0.866025403784439  "," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[1]),str(unpaired[2])," b ",'  '.join(map(str, paired)),str(unpaired[3]), file=flist)
    return 3

def ne5spin2(paired,unpaired,irun):
    flist = open('list'+str(irun)+'.txt','a')
    npaired = len(paired)
    nunpaired = len(unpaired)
    if(not(nunpaired==5)):
        print("error in ne5spin1, I stop")
        sys.exit()

    print(4, file=flist)
    print(0.5," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[2]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[1]), str(unpaired[3]), file=flist)
    print(0.5," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[3]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[1]), str(unpaired[2]), file=flist)
    print(0.5," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[2]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[3]), file=flist)
    print(0.5," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[3]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[2]), file=flist)
    print(6, file=flist)
    print("0.57735026918962584"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[1]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[2]), str(unpaired[3]), file=flist)
    print("0.57735026918962584"," a ",'  '.join(map(str, paired)),str(unpaired[2]),str(unpaired[3]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[1]), file=flist)
    print("0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[2]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[1]), str(unpaired[3]), file=flist)
    print("0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[3]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[2]), file=flist)
    print("-0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[3]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[1]), str(unpaired[2]), file=flist)
    print("-0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[2]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[3]), file=flist)
    print(6, file=flist)
    print("0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[3]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[1]), str(unpaired[2]), file=flist)
    print("0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[2]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[1]), str(unpaired[3]), file=flist)
    print("0.57735026918962584"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[1]),str(unpaired[2])," b ",'  '.join(map(str, paired)), str(unpaired[1]), str(unpaired[4]), file=flist)
    print("-0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[3]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[2]), file=flist)
    print("-0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[2]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[3]), file=flist)
    print("-0.57735026918962584"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[1]),str(unpaired[2])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[4]), file=flist)
    print(6, file=flist)
    print("0.40824829046386302"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[1]),str(unpaired[2])," b ",'  '.join(map(str, paired)), str(unpaired[3]), str(unpaired[4]), file=flist)
    print("0.40824829046386302"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[1]),str(unpaired[3])," b ",'  '.join(map(str, paired)), str(unpaired[2]), str(unpaired[4]), file=flist)
    print("0.40824829046386302"," a ",'  '.join(map(str, paired)),str(unpaired[2]),str(unpaired[3]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[1]), file=flist)
    print("-0.40824829046386302"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[2]),str(unpaired[3])," b ",'  '.join(map(str, paired)), str(unpaired[1]), str(unpaired[4]), file=flist)
    print("-0.40824829046386302"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[2]),str(unpaired[3])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[4]), file=flist)
    print("-0.40824829046386302"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[1]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[2]), str(unpaired[3]), file=flist)
    print(6, file=flist)
    print("0.57735026918962584"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[1]),str(unpaired[3])," b ",'  '.join(map(str, paired)), str(unpaired[2]), str(unpaired[4]), file=flist)
    print("-0.57735026918962584"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[1]),str(unpaired[2])," b ",'  '.join(map(str, paired)), str(unpaired[3]), str(unpaired[4]), file=flist)
    print("0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[2]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[1]), str(unpaired[3]), file=flist)
    print("0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[2]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[3]), file=flist)
    print("-0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[0]),str(unpaired[3]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[1]), str(unpaired[2]), file=flist)
    print("-0.28867513459481292"," a ",'  '.join(map(str, paired)),str(unpaired[1]),str(unpaired[3]),str(unpaired[4])," b ",'  '.join(map(str, paired)), str(unpaired[0]), str(unpaired[2]), file=flist)
    return 5

