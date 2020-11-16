#!/usr/bin/env python
import sys
import argparse

import numpy as np
import tensorflow as tf
from tensorflow import keras

aa =['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H',
 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class Seq2FeatVec:
    def __init__(self,fname, dname, train_type=None):
        self.infile = fname
        self.wrkdir = dname
        self.train_type=train_type
        if train_type==None:
            self.train_type = dname
                 
        self.Dexons = {
            "igfish":['back','exon1_M','exon2_M','exon3_M','exon4_M','exon1_D',
                       'exon1B_D','exon2_D','exon3_D','exon4_D','exon5_D','exon6_D',
                       'exon2_T','exon3_T','exon4_T'],
            "igreptiles": ['back','exon1_M','exon2_M','exon3_M','exon4_M', 'exon1_D',
                            'exon2_D','exon3_D','exon4_D','exon5_D','exon6_D','exon7_D',
                            'exon8_D','exon9_D','exon10_D','exon11_D','exon1_Y','exon2_Y',
                            'exon3_Y','exon4_Y','exon1_A','exon2_A','exon3_A','exon1_E',
                            'exon2_E','exon3_E','exon4_E','exon1_XA','exon2_XA','exon3_XA',
                            'exon4_XA'], 
            "mhc":['back','exon1_I','exon2_I','exon3_I','exon2_DA','exon3_DA',
                   'exon2_DB','exon3_DB'],
            "vsfish":['back','ighv','igkv','iglv','trav','trbv','trgv','trdv','nitr'],
            "chsfish":['back','exon1_M','exon2_M','exon3_M','exon4_M','exon1_D','exon1B_D','exon2_D',
                       'exon3_D','exon4_D','exon5_D','exon6_D','exon1_T','exon2_T','exon3_T','exon4_T']
        }

        
    def seq2feature_vec(self):
        X0 = []
        Y0 =[]
        for w in range(len(self.Dexons[self.train_type])):
            count = 0
            for record in SeqIO.parse(self.wrkdir + "/" + self.infile,'fasta'):
                gro = record.id.split('-')[0]

                if "Seriola" in record.id:
                    print(record.id.split('-'))
                #if (gro==self.Dexons[self.train_type][w]):
                #    print(len(record.seq))
                if self.Dexons[self.train_type][w] == gro and len(record.seq)>=42:
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
                    Y0.append(w)
                    count +=1

            print(self.Dexons[self.train_type][w],'-----> ',count)


        X = np.asarray(X0)
        print(X.shape)
        Y = np.asarray(Y0)
        print(Y.shape)

        P =[]
        PP =[]
        co = 0
        for x in X:
            co +=1
            if co%3 == 0:
                PP.append(x)
            else:
                P.append(x)
        Q =[]
        QQ =[]
        co = 0
        for y in Y:
            co +=1
            if co%3 == 0:
                QQ.append(y)
            else:
                Q.append(y)            

        x_train = np.array(P,dtype=np.int)
        x_test  = np.array(PP,dtype=np.int)
        y_train = np.array(Q,dtype=np.int)
        y_test  = np.array(QQ,dtype=np.int)

        np.save(self.wrkdir + "/x_train",x_train)
        np.save(self.wrkdir + "/y_train",y_train)
        np.save(self.wrkdir + "/x_test",x_test)
        np.save(self.wrkdir + "/y_test",y_test)


    def ml_accuracy(self):
        x_train=np.load(self.wrkdir + "/x_train.npy")
        y_train=np.load(self.wrkdir + "/y_train.npy")
        x_test=np.load(self.wrkdir + "/x_test.npy")
        y_test=np.load(self.wrkdir + "/y_test.npy")        
        model = keras.Sequential([
            keras.layers.Dense(128, activation=tf.nn.relu),
            keras.layers.Dense(len(self.Dexons[self.train_type]), activation=tf.nn.softmax)
        ])
        model.compile(optimizer='adam',
                      loss='sparse_categorical_crossentropy',
                      metrics=['accuracy'])
        model.fit(x_train, y_train, epochs=10)
        est_loss, test_acc = model.evaluate(x_train, y_train)
        print('Test accuracy:', test_acc)


        
    def run(self):
        print("Feat for training:")
        print(" ... Exons=", self.Dexons[self.train_type])
        print(" ... infile=",self.infile)
        self.seq2feature_vec()
        self.ml_accuracy()

        
def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Training Code')
    parser.add_argument('-i', '--infasta',
                        help='The input fasta file',
                        required=True,
                        default='input.fasta')                            
    parser.add_argument('-d', '--dirname',
                        help='The working directory',
                        required=True,
                        default="")
    parser.add_argument('-t', '--traintype',
                        help='Training Type',
                        required=False,
                        default="igfish")
    results = parser.parse_args(args)
    return (results.infasta, results.dirname, results.traintype)

# -----------------------------------------------
if __name__ == '__main__':
    infasta, dirname, traintype=  check_arg(sys.argv[1:])
    S = Seq2FeatVec(infasta, dirname)
    S.run()

    
