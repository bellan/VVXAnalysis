######################################################
# A python helper module, written for VZZAnalyzer,   #
# to unpickle (load) scikit-learn objects            #
#                                                    #
# Author: A. Mecca (alberto.mecca@edu.unito.it)      #
######################################################

from os.path import exists
from os import environ
from sklearn.ensemble import AdaBoostClassifier
from pickle import Unpickler
from warnings import filterwarnings

def test():
	print('exists("src")', exists("src"))


def load_object(path):
	if(not exists(path)):
		print('Error: path "%s" does not exist.' % (path))
		return None
	with open(path, "rb") as file:
		unpickler = Unpickler(file)
		filterwarnings("ignore", category=UserWarning) # suppress version mismatch warnings
		loaded_Ada_VZZ = unpickler.load()
		filterwarnings("default") # restore default behaviour
		del unpickler
	return loaded_Ada_VZZ
	
def predict(ADA, list_data):
	#print("Predicting. type(list_data):", type(list_data))
	#print("Prediction: ", ADA.predict_proba([list_data]))
	res = ADA.predict_proba([list_data])[:, 1][0]
	#print("type(res) =", type(res))
	#print("Returning a %s  (value = %.4f)" % (type(res), res))
	return res;
	
def pyfloat_from_prediction(pred):
	return float(pred[0, 1])

def test_ADA(ADA):
	print("is ADA instance of AdaBoostClassifier?", isinstance(ADA, AdaBoostClassifier))

def test_list(list):
	print("list:", list)
	
def print_type(obj):
	print(type(obj))
	
def print_obj(obj):
	print(obj)

if(__name__ == "__main__"):
	ll = load_object()
	test(ll)
	X = [0.1, -3.5, 12.6, 32., 0.5]
	pred = predict( ll, [X] )
	print("prediction = 0: %.1f  1: %.1f" % tuple(pred[0]) )

