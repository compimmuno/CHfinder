#!/usr/bin/env python

import argparse
import itertools
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re
import pickle 
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIWWW
from time import time
import sys, os
import tensorflow as tf
from tensorflow import keras
import numpy as np
import matplotlib.pyplot as plt

aa =['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L',
     'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

class ExonFinder:
    def __init__(self, spfile, spdir, traindir, train_type=None, query_file=None):
        self.spfile = spfile
        self.spdir = spdir
        self.traindir = traindir
        self.train_type=train_type
        self.query_file = query_file
           
        self.Dexons = {
            "igfish":['back','exon1_M','exon2_M','exon3_M','exon4_M','exon1_D',
                       'exon1B_D','exon2_D','exon3_D','exon4_D','exon5_D','exon6_D',
                       'exon2_T','exon3_T','exon4_T'],
            "igreptiles": ['back','exon1_M','exon2_M','exon3_M','exon4_M', 'exon1_D',
                            'exon2_D','exon3_D','exon4_D','exon5_D','exon6_D','exon7_D',
                            'exon8_D','exon9_D','exon1_D','exon11_D','exon1_Y','exon2_Y',
                            'exon3_Y','exon4_Y','exon1_A','exon2_A','exon3_A','exon1_E',
                            'exon2_E','exon3_E','exon4_E','exon1_XA','exon2_XA','exon3_XA',
                            'exon4_XA'], 
            "mhc":['back','exon1_I','exon2_I','exon3_I','exon2_DA','exon3_DA',
                   'exon2_DB','exon3_DB'],
            "chsfish":['back','exon1_M','exon2_M','exon3_M','exon4_M','exon1_D','exon1B_D','exon2_D',
                       'exon3_D','exon4_D','exon5_D','exon6_D','exon1_T','exon2_T','exon3_T','exon4_T']
        }


    def run_tblastn(self):
        #tbin = "/home/david/anaconda3/envs/bioinf/bin/tblastn"
        tbin = "tblastn"
        outfile = self.spdir + "/tabla.txt"
        species_file = self.spdir+"/"+self.spfile
        queryfile= self.traindir + "/" + self.query_file
        blastn_cline = NcbiblastxCommandline(tbin,
                                             query = queryfile, 
                                             subject = species_file, 
                                             evalue = 0.001,
                                             out = outfile, outfmt =6)
        print(blastn_cline)
        stdout, stderr = blastn_cline()


    def process_tblastn_table(self):
        tablafile = self.spdir+'/tabla.txt'
        A = open(tablafile,'r')
        ide = []
        contigs = []
        for line in A:
            par = line.split('\t')
            if '#' not in line:
                ID = par[1]
                a = par[8]
                b = par[9]
                if int(b) - int(a) >=0:
                    frame = 1
                    start = a
                else:
                    frame = -1
                    start = b
                hit = ID,frame,start
                contigs.append(hit)
                if ID not in ide:
                    ide.append(ID)
        return contigs, ide

                    
        
    def candidate_exons(self, contigs, ide):
        ofile = open(self.spdir+'/aexones.fasta', "w")
        oofile = open(self.spdir+'/nexones.fasta', "w")

        species_file = self.spdir+"/"+self.spfile        
        
        for record in SeqIO.parse(species_file, "fasta"):
            ee = record.id
            if ee in ide:
                seq = record.seq.upper()
                print (ee,'long seq ',len(seq)) 
                ini1 = [i.start() for i in re.finditer('AG',str(seq))]

                ini=[]
                for m,n,p in contigs:
                    if int(n) >= 0 and m == ee:
                        w = int(p) -4000
                        ww = int(p) + 4000
                        for b in ini1:
                            if b >=w and b<= ww:
                                if b not in ini:
                                    ini.append(b)



                fin1 = [i.start() for i in re.finditer('GT',str(seq))]
                fin = []
                for m,n,p in contigs:
                    if int(n) >= 0 and m == ee:
                        w = int(p) -4000
                        ww = int(p) + 4000
                        for b in fin1:
                            if b >=w and b<= ww:
                                if b not in fin:
                                    fin.append(b)
                                    
                se = [(i,j) for i,j in itertools.product(ini,fin) if j>i+240 and j<i+390]
                print ('secuencias posibles foward',len(se)) 
                for x in range(len(se)):
                    a = se[x][0]
                    b = se[x][1]
                    exon = seq[a+2:b]       
                    tras = exon[2:]
                    p = tras.translate(to_stop=True)
                    if len(p)>=80 :   #and not '*' in p and not 'X' in p:
                        nrec = SeqRecord(exon, id = ee + '-' + str(a) + '-' + str(b) ,description = 'plus')
                        arec = SeqRecord(p, id = ee + '-' + str(a) + '-' + str(b)  ,description = 'plus')
                        ofile.write(arec.format('fasta'))
                        oofile.write(nrec.format('fasta'))


                seq = record.seq.reverse_complement()
                seq = seq.upper()
                ancho = len(seq)
                ini1 = [i.start() for i in re.finditer('AG',str(seq))]

                ini=[]
                ini=[]
                for m,n,p in contigs:
                    if int(n) <= 0 and m == ee:
                        pp = ancho-int(p)
                        w = pp -4000
                        ww = pp + 4000
                        for b in ini1:
                            if b >=w and b<= ww:
                                if b not in ini:
                                    ini.append(b)

                fin1 = [i.start() for i in re.finditer('GT',str(seq))]
                fin = []
                for m,n,p in contigs:
                    if int(n) <= 0 and m == ee:
                        pp = ancho-int(p)
                        w = pp -4000
                        ww = pp + 4000
                        for b in fin1:
                            if b >=w and b<= ww:
                                if b not in fin:
                                    fin.append(b)

                se = [(i,j) for i,j in itertools.product(ini,fin) if j>i+230 and j<i+390 and ((np.abs(j-i)-2)%3==0)]
                print ('secuencias posibles reverse',len(se)) 
                for x in range(len(se)):
                    a = se[x][0]
                    b = se[x][1]              
                    exon = seq[a+2:b]            
                    tras = exon[2:-1]
                    p = tras.translate(to_stop=True)
                
                    if len(p)>=80 and not '*' in p and not 'X' in p:
                        aa = ancho-b
                        bb = ancho-a
                        nrec = SeqRecord(exon, id = ee + '-' + str(aa) + '-' + str(bb) ,description = 'minus')
                        arec = SeqRecord(p, id = ee + '-' + str(aa) + '-' + str(bb),description = 'minus')
                        ofile.write(arec.format('fasta'))
                        oofile.write(nrec.format('fasta'))


        ofile.close() 
        oofile.close() 



    ### (latest variation of code for CHfinder).
    def candidate_exons_chfinder(self, contigs, ide):
        ofile = open(self.spdir+'/aexones.fasta', "w")
        oofile = open(self.spdir+'/nexones.fasta', "w")

        species_file = self.spdir+"/"+self.spfile        
        #animal="ANIMAL"
        animal=os.path.basename(self.spdir)

        
        rep = []
        for record in SeqIO.parse(species_file, "fasta"):
            ee = record.id
            f = 0
            r = 0
            if ee in ide:
                seq = record.seq
                seq = seq.upper()
                for m, n, p in contigs:
                    if int(n) >= 0 and m == ee:
                        yy = int(p) - 100
                        k = seq[yy:yy + 1000]
                        ini1 = [i.start() for i in re.finditer('AG', str(k))]
                        fin1 = [i.start() for i in re.finditer('GT', str(k))]
                        se = [(i, j) for i, j in itertools.product(ini1, fin1) if
                              j > i + 240 and j < i + 390 and ((np.abs(j - i) - 2) % 3 == 0)]
                        DF = []
                        for r, s in se:
                            x = r + yy
                            y = s + yy
                            pis = x, y
                            if pis not in DF:
                                DF.append(pis)
                                exon = seq[x + 2:y]
                                tras = exon[2:]
                                p = tras.translate(to_stop=True)
                                if len(p) >= 80:
                                    nrec = SeqRecord(exon, id=animal + '-' + ee + '-' + str(x + 4) + '-' + str(y),
                                                     description='plus')
                                    arec = SeqRecord(p, id=animal + '-' + ee + '-' + str(x + 4) + '-' + str(y),
                                                     description='plus')
                                    if arec.id not in rep:
                                        ofile.write(arec.format('fasta'))
                                        oofile.write(nrec.format('fasta'))
                                        rep.append(arec.id)
                                        f += 1
                print(animal, '----> ', record.id, len(record.seq))
                print('foward squences posibles ', f)
                seq = record.seq.reverse_complement()
                seq = seq.upper()
                ancho = len(seq)

                for m, n, p in contigs:
                    if int(n) <= 0 and m == ee:
                        pp = ancho - int(p)
                        yy = int(pp) - 1000
                        k = seq[yy:yy + 2000]
                        ini1 = [i.start() for i in re.finditer('AG', str(k))]
                        fin1 = [i.start() for i in re.finditer('GT', str(k))]
                        se = [(i, j) for i, j in itertools.product(ini1, fin1) if
                              j > i + 240 and j < i + 390 and ((np.abs(j - i) - 2) % 3 == 0)]
                        DR = []
                        for r, s in se:
                            x = r + yy
                            y = s + yy
                            pis = x, y
                            if pis not in DR:
                                DR.append(pis)
                                exon = seq[x + 2:y]
                                tras = exon[2:]
                                p = tras.translate(to_stop=True)
                                if len(p) >= 80:
                                    aa = ancho - y
                                    bb = ancho - x
                                    nrec = SeqRecord(exon, id=animal + '-' + ee + '-' + str(aa) + '-' + str(bb),
                                                     description='minus')
                                    arec = SeqRecord(p, id=animal + '-' + ee + '-' + str(aa) + '-' + str(bb),
                                                     description='minus')
                                    if arec.id not in rep:
                                        ofile.write(arec.format('fasta'))
                                        oofile.write(nrec.format('fasta'))
                                        rep.append(arec.id)
                                        r += 1

                print('reverse squences posibles ', r)

        ofile.close()
        oofile.close()



    ### --------------------


    def predict_viable_exons(self):        
        x_train=np.load(self.traindir + "/x_train.npy")
        y_train=np.load(self.traindir + "/y_train.npy")
        x_test=np.load(self.traindir + "/x_test.npy")
        y_test=np.load(self.traindir + "/y_test.npy")        

        model = keras.models.Sequential([
            keras.layers.Dense(128, activation='relu'),
            tf.keras.layers.Dropout(0.2),
            keras.layers.Dense(len(self.Dexons[self.train_type]), activation='softmax')
        ])
        model.compile(optimizer='adam',
                      loss='sparse_categorical_crossentropy',
                      metrics=['accuracy'])
        model.fit(x_train, y_train, epochs=5)
        test_loss, test_acc = model.evaluate(x_train, y_train)
        print('Test accuracy:', test_acc)

        bas = self.spdir+'/aexones.fasta'

        
        ID = []
        X0 = []
        contigs = []
        count = 0
        for record in SeqIO.parse(bas, "fasta"):
            ID.append(record.description)
            par = record.id.split('-')
            if par[0] not in contigs:
                contigs.append(par[0])
            seq = record.seq[:40]+record.seq[-40:]
            pp = record.seq
            fre = []
            for x in range(len(seq)):
                a = seq[x]
                for y in aa:
                    if a == y:
                        fre.append(1)
                    else:
                        fre.append(0)
            for x in range(len(aa)):
                for y in range(len(aa)):
                    i = aa[x] + aa[y]
                    q = pp.count(i)
                    fre.append(q)
            X0.append(fre)
            count +=1

        print('secuencias  ',count)    
        X = np.asarray(X0)
        print(X.shape)

        X0 = np.array(X,dtype=np.int)
        print(X0.shape)
        predictions = model.predict(X0)
        print (predictions[0])
        print('numerical transform')
        oo = open(self.spdir+'/aexones_asignados.txt','w')
        for a in range(len(predictions)):
            pp = round(max(predictions[a]),3)
            p = ID[a]+'-----------> '+str(pp)+'-'+self.Dexons[self.train_type][np.argmax(predictions[a])]+'\n'
            oo.write(p)
        oo.close()
        print('asignacion hecha')

        good = []
        count = 0
        nlp = []
        for y in range(len(predictions)):
            locus = self.Dexons[self.train_type][np.argmax(predictions[y])]
            if locus != 'back':
                pro = np.max(predictions[y])
                prob = round(pro,3)
                name = ID[y]
                p = name,locus,prob
                nlp.append(p)

        for b in contigs:
            exonesplus = []
            exonesminus = []

            for x,y,z in nlp:
                st = x.split('-')
                p = x,y,z,st[2]    ### changed this....
                print(p)
                
                if b in x and y != 'back' and 'plus' in x:
                    exonesplus.append(p)
                if b in x and y != 'back' and 'minus' in x:
                    exonesminus.append(p)

            exp = sorted(exonesplus,key=lambda x: x[3])
            exm = sorted(exonesminus,key=lambda x: x[3]) 

            #print(exp)
            #k = range(200)
            if len(exp)>=1:
                #print(exp[0].split("-"))
                k = range(int(exp[0][3])-2,int(exp[0][3])+400)
                partes =[]
                for z in range(len(exp)):
                    if int(exp[z][3]) in k:
                        partes.append(exp[z])
                    elif z == len(exp):
                        zz = max(partes, key=lambda t: t[2])
                        if zz not in good:
                            good.append(zz)
                            count +=1
                    else:
                        if len(partes)>= 1:
                            zz = max(partes, key=lambda t: t[2])
                            if zz not in good:
                                good.append(zz)
                                count +=1
                        if len(exp)>= 1:
                            k = range(int(exp[z][3])-2,int(exp[z][3])+400)
                            partes = []
                            partes.append(exp[z])

                if len(partes)>= 1:
                    zz = max(partes, key=lambda t: t[2])
                    if zz not in good:
                        good.append(zz)
                        count +=1


            #k = range(200)
            if len(exm)>=1:
                k = range(int(exm[0][3])-2,int(exm[0][3])+400)
                partes =[]
                for z in range(len(exm)):
                    if int(exm[z][3]) in k:
                        partes.append(exm[z])
                    else:  
                        if len(partes)>= 1:
                            zz = max(partes, key=lambda t: t[2])
                            if zz not in good:
                                good.append(zz)
                                count +=1
                        if len(exm)>= 1:
                            k = range(int(exm[z][3])-2,int(exm[z][3])+400)
                            partes = [exm[z]]
                #print(partes)            
                if len(partes)>= 1:
                    zz = max(partes, key=lambda t: t[2])
                    if zz not in good:
                        good.append(zz)
                        count +=1
                        
        print ('good ' ,count)  

        ccc = 0 
        o = open(self.spdir+'/result-'+self.train_type+'.fasta','w')       
        for re in SeqIO.parse(bas,'fasta'):
            for a,b,c,d in good:
                pir = a.split(' ')
                if re.id in pir[0] and float(c) >= 0.2:
                    rec = SeqRecord(re.seq,id = pir[0]+'-'+b+'-'+str(c),description = pir[1])


                    o.write(rec.format('fasta'))
                    ccc +=1
        o.close()                    
        print ('numero de secuencias exones ',ccc)


        
    # ---------
    def run(self):
        print("Input parameters for prediction:")
        print(" ... Exons=", self.Dexons[self.train_type])
        print(" ... speciesfile=",self.spfile)
        print(" ... speciesdir=",self.spdir)
        print(" ... traindir=",self.traindir)
        print(" ... traindir=",self.query_file)        


        # step 1.  Run tbalstn to identify location of hits 
        self.run_tblastn()

        # step 2. process hit file
        contigs, ide = self.process_tblastn_table()
        print(contigs)



        # step 3. select candidate sequences from genome
        #self.candidate_exons(contigs, ide)
        self.candidate_exons_chfinder(contigs, ide)


        # step 4: run prediction
        self.predict_viable_exons()




#------------------------------
def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Exon Finders')
    parser.add_argument('-i', '--speciesfile',
                        help='The species fasta file',
                        required=True,
                        default='input.fasta')                            
    parser.add_argument('-d', '--speciesdir',
                        help='The species directory',
                        required=True,
                        default="")
    parser.add_argument('-T', '--trainingdir',
                        help='The training directory',
                        required=True,
                        default="")
    parser.add_argument('-E', '--exontype',
                        help='Training Type',
                        required=False,
                        default="igfish")
    parser.add_argument('-q', '--queryfile',
                        help='Query File',
                        required=False,
                        default="qfile.fasta")


    results = parser.parse_args(args)
    return (results.speciesfile, results.speciesdir, results.trainingdir, results.exontype, results.queryfile)

# -----------------------------------------------
if __name__ == '__main__':
    speciesfile, speciesdir, traindir, exontype, queryfile =  check_arg(sys.argv[1:])
    #query = 'query-reptiles.fasta.txt'


    S = ExonFinder(speciesfile, speciesdir, traindir, exontype, queryfile)
    S.run()
