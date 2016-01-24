#!/usr/bin/python
# -*- coding: utf-8 -*-
import time
import numpy as np
from math import log, ceil
import bitarray as ba
from ast import literal_eval as l_eval


##### FUNCIONS #####

##########################################
##               HUFFMAN                ##
##########################################
def huffman(ifile):
    if ifile[len(ifile)-3:len(ifile)]!='huf':
        text=open(ifile,'r').read()
        # Frequency of each char
        cf_dic = {}
        for char in text:
            cf_dic[char] = cf_dic.get(char, 0) + 1

        # Convert cf_dict to list of tuples, flip elements and order highest freq first
        value_key = sorted([(f, k) for k, f in cf_dic.items()], reverse=True)

        # Convert frequency to probability
        freq_total= float(sum([freq[0] for freq in value_key]))
        value_key=zip([freq[0]/freq_total for freq in value_key],[char[1] for char in value_key])

        bolean=True
        tmp=value_key
        combinacions=[]

        # Huffman algorithm
        while bolean:
            combinacions.append([tmp[-2][1],'0'])
            combinacions.append([tmp[-1][1],'1'])
            tmp2=tmp
            tmp=tmp[:-2]
            tmp.append((tmp2[-1][0]+tmp2[-2][0],tmp2[-1][1]+tmp2[-2][1]))
            tmp.sort(key=lambda tup:tup[0], reverse=True)
            if len(tmp)==1:
                bolean=False

        taula={}  # Taula a passar
        for comb,b in combinacions:
            for c in list(comb):
                try:
                    taula[c]= b+ taula[c]
                except:
                    taula[c]=b

        textc=[]
        taulaint={k: map(int,v) for k, v in taula.items()}
        for char in text:
            textc+=taulaint[char]
        
        textb=ba.bitarray(textc)
        #Taula en bitarray
        tableb=ba.bitarray()
        tableb.frombytes(str(taula))

        # Primers bytes que ens diran la longitud de la taula i els bits que he afegit al final (que no pertanyen a la compressió, sino a acabar el byte)
        lentableb=bin(len(tableb))[2:]
        aux=ba.bitarray('0'*(16-len(lentableb))+lentableb)
        bitsafegits=len(textb) % 8
        aux+=ba.bitarray('0'*(8-len(bin(bitsafegits)[2:]))+bin(bitsafegits)[2:])
        tableb=aux+tableb
        # Write file 
        out=open(ifile+'.huf','wb')
        tableb.tofile(out)
        textb.tofile(out)
    # Extract
    else:
        data=np.unpackbits(np.fromfile(ifile,dtype='uint8')).tolist()
        lendic=int(''.join(map(str,data[:16])),2)
        bitsafegits=int(''.join(map(str,data[16:24])),2)
        
        tableb=ba.bitarray(data[24:(24+lendic)])
        table=l_eval(tableb.tobytes())
        #Invert the table
        table={ v: k for k, v in table.items()}
        data=data[(24+lendic):]
        text=''
        bits=''
        for b in data:
            bits+=str(b)
            if bits in table:
                text+= table[bits]
                bits=''

        out=open('HUF_'+ifile[:-4],'w')
        out.write(text)


##########################################
##             LEMPEL-ZIV               ##
##########################################
def lz(ifile):
    if ifile[-3:]!='zip':
        # Comprimir arxiu
        data=np.unpackbits(np.fromfile(ifile,dtype='uint8')).tolist()
        index=1
        order=[]
        bits=1
        lastbits=1
        dic={1:[0,0,1]} # Trick to avoid errors with 0s in front of 1s. It will have to be filtered
        for b in data:
            bits=bits << 1 | b
            if bits not in dic:
                dic[bits]=[index,dic[lastbits][0],b]
                index+=1
                order.append(bits)
                bits=1
                lastbits=1
            lastbits=bits

        if bits!=1: # Per incloure els bits que ens puguin quedar orfes
            order.append(bits)
                
        #Número de bits per cada pointer+últim bit
        nbits=int(ceil(log(len(dic)-1)/log(2.)))+1 

        # Els primers 8 bits ens indicaran la longitud de cada codi Nbits+1 i els segons 8 bits ens indicaran el numero de bits afegits al final de més (en cas que 
        
        textb=ba.bitarray()
        textb=[]
        for phr in order:
            bit=dic[phr]
            textb+=[0]*(nbits-len(bin(bit[1])[2:])-1)+map(int,bin(bit[1])[2:])+[bit[2]]
        textb=ba.bitarray(textb)
        textc=ba.bitarray('0'*(8-len(bin(nbits)[2:]))+bin(nbits)[2:])
        bitsafegits=len(textb) % 8
        textc+=ba.bitarray('0'*(8-len(bin(bitsafegits)[2:]))+bin(bitsafegits)[2:])

        textc+=textb
        out=open(ifile+'.zip','wb')
        textc.tofile(out)

    
    # Descomprimir arxiu
    else:

        data=np.unpackbits(np.fromfile(ifile,dtype='uint8')).tolist()
        dic={}
        bits=None
        lenbits=int(''.join(map(str,data[:8])),2)
        bitsafegits=int(''.join(map(str,data[8:16])),2)
        # Packing the list in lists of lenbits elements
        data=data[16:]
        data=[data[i:i+lenbits] for i in range(0,len(data)) if (i%lenbits)==0]
        j=1
        for bits in data[:-1]: 
            indexant=int(''.join(map(str,bits[:-1])),2)
            dic[j]=[indexant,bits[-1]]
            j+=1

        text=[]
        t1=time.clock()
        # Schema: index:phrase
        phrases={}
        for i in xrange(1,len(dic)+1):
            try:
                phrase=phrases[dic[i][0]]+[dic[i][1]]
                phrases[i]=phrase
            except:
                phrases[i]=[dic[i][1]]
            text+=phrases[i]
        t2=time.clock()

        textb=ba.bitarray(text)
        ofile=ifile

        out=open('LZ_'+ifile[:-4],'wb')
        textb.tofile(out)

##########################################
##             LEMPEL-ZIV               ##
##########################################
def equals(orig,test):
    o=open(orig,'r').read()
    t=open(test,'r').read()
    return o==t


##### PROGRAMA #####

print 'Choose what do you want to do:'
print '\t h \t Compress/extract with huffman algorithm'
print '\t l \t Compress/extract with LX78 algorithm'
print '\t c \t Check if the extracted files are equal'
print '\t t \t Together'
action=raw_input()
print 'Write the input filename:'
ifile=raw_input()


if action=='t':
    try:
        t1=time.clock()
        huffman(ifile)
        t2=time.clock()
        huffman(ifile+'.huf')
        t3=time.clock()
        lz(ifile)
        t4=time.clock()
        lz(ifile+'.zip')
        t5=time.clock()
        print 'HUFFMAN'
        print '\tCompression:',t2-t1
        print '\tExtraction:',t3-t2
        print '\tTotal:', t3-t1
        print '\tEqual?', equals(ifile,'HUF_'+ifile)
        print 
        print 'LEMPEL-ZIV'
        print '\tCompression:', t4-t3
        print '\tExtraction:', t5-t4
        print '\tTotal:', t5-t3
        print '\tEqual?', equals(ifile,'LZ_'+ifile)
    except:
        print '\n\n ERROR'
    

else:
    try:
        t1=time.clock()
        if action=='h': 
            huffman(ifile)
        elif action=='l':
            lz(ifile)
        elif action=='c':
            print 'Huffman:', equals(ifile,'HUF_'+ifile)
            print 'LZ:', equals(ifile,'LZ_'+ifile)
        else:
            print 'I did not underestand it. Write only h or l'
        t2=time.clock()
        print 'Time elapsed: ', t2-t1

    except:
        print '\n Error'


