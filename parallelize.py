#!/usr/bin/env python

import subprocess
import os
import threading
import Queue
import time

class Worker(threading.Thread):
	def __init__(self, queue):
		self.__queue = queue
		threading.Thread.__init__(self)
	
	def run(self):
		while True:
			item = self.__queue.get()
			if item is None:
				break
			#w = subprocess.Popen("./waitAndPrint.py "+str(item)+"Hi", shell=True)
			#w.wait()
			time.sleep(.1)
			print "Done %s"%item

lock = threading.Lock()
class LockedWorker(threading.Thread):
	def __init__(self, queue):
		self.__queue = queue
		threading.Thread.__init__(self)
	
	def run(self):
		while True:
			lock.acquire()
			item = self.__queue.get()
			if item is None:
				lock.release()
				break
			while os.path.exists("globallock.txt"):
				print "Waiting for global lock to be released"
				time.sleep(0.5)
			open("globallock.txt", "w")
						
			w = subprocess.Popen("./updateCSV.py "+str(item)+"Hi", shell=True)
			w.wait()
			print "Done %s"%item
			os.remove("globallock.txt")
			lock.release()

NUMWORKERS = 2
text = "1 2 3 4 5 6 7 8 9 10".split(' ')
NUMJOBS = len(text)

if __name__ == '__main__' :
	queue = Queue.Queue(4)
	for i in xrange(NUMWORKERS):
		LockedWorker(queue).start()
	
	for item in text:
		queue.put(item)
	
	for item in xrange(NUMWORKERS):
		queue.put(None)

else:
	print "Not main, not running"
