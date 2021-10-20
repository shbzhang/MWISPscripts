import glob, re
import numpy as np
from math import floor,ceil

class CellMap(np.ndarray):
	'''
	Check the number of files in MWISP cells in terminal
	Parameters
	----------
	patterns : str or list
		a pathname with '*' or '?' to match files
	l,b : list, optional
		2-elements list that define the map range to show, default shows cells with l=[0,10],b=[-5,5]
	deluxe : bool, optional
		format of output, set this keyword true to show a deluxe format.
	reg : str or re.Pattern, optional
		format of cell name.
	Versions
	--------
	Oct,20,2021,v1.0
		a python version of mapfile.pro, using a subclass of numpy.ndarray
	Examples
	--------
	#check observed cells
	>>>havecube = CellMap('/share/data/mwisp/G???+00/*.fits',l=[0,20])
	>>>havecube

	#change position
	>>>havecube.l = [35,50]
	>>>havecube.b = [-6,6]
	>>>havecube.deluxe = True
	>>>havecube
	
	#search of missing info file
	>>>haveinfo = CellMap('/share/data/mwisp/infofiles/*info.txt')
	>>>missinginfo = (havecube>0) & (haveinfo==0)
	>>>missinginfo
	>>>missinginfo.where(all=False)
	'''
	def __new__(cls, *patterns, l=[0,10], b=[-5,5], \
		deluxe=False, reg=re.compile('([012]\d\d[05])([-\+][01]\d[05])')):
		obj = np.asarray(np.zeros((81,360*2),dtype=np.int32)).view(cls)
		obj.reg = reg
		obj.patterns = patterns
		obj.l = l
		obj.b = b
		obj.deluxe = deluxe
		return obj

	def __array_finalize__(self, obj):
		if obj is None: return
		for key in ['reg','patterns','l','b','deluxe']:
			self.__dict__[key] = getattr(obj, key, None)

	def __setattr__(self, key, value):
		self.__dict__[key] = value
		if (key == 'patterns') and (value is not None):
			if type(value) is str: self.__dict__[key] = (value,)
			self._count()

	@property
	def laxis(self):
		return np.arange(3595,-1,-5)

	@property
	def baxis(self):
		return np.arange(200,-201,-5)

	@property
	def lindex(self):
		l0 = (self.laxis[0] - floor(max(self.l)*2)*5) //5
		l1 = (self.laxis[0] - ceil(min(self.l)*2)*5) //5
		return (l0,l1+1)

	@property
	def bindex(self):
		b0 = (self.baxis[0] - floor(max(self.b)*2)*5) //5
		b1 = (self.baxis[0] - ceil(min(self.b)*2)*5) //5
		return (b0,b1+1)

	def _count(self):
		#search files
		files = []
		for pt in self.patterns:
			files += glob.glob(pt)
		#find cellname in files
		self[:] = 0
		for file in files:
			found = re.search(self.reg, file)
			if found is not None:
				lname,bname = found.group(1,2)
				lidx = (self.laxis[0]-int(lname))//5
				bidx = (self.baxis[0]-int(bname))//5
				self[bidx,lidx] += 1

	def __str__(self):
		lindex, bindex = self.lindex, self.bindex
		laxis = self.laxis[lindex[0]:lindex[1]]/10
		baxis = self.baxis[bindex[0]:bindex[1]]/10
		count = self[bindex[0]:bindex[1], lindex[0]:lindex[1]]
		if self.deluxe:
			#head
			hline = '-----' + '+---'*count.shape[1] + '+-----\n'
			output = r'  b\l'
			if laxis[0]%1==0:
				for l in laxis[0::2]: output += '|%-7.1f' % l
			else:
				output += ' '*4
				for l in laxis[1::2]: output += '|%-7.1f' % l
			output += '\n'
			output += hline
			#body
			for b,line in zip(baxis, count):
				output += '%+5.1f' % b
				for value in line: output += '|%2i ' % value if value!=0 else '|   '
				output += '|%+-5.1f\n' % b
				output += hline
			#rear
			output += ' '*5
			if laxis[0]%1==0:
				output += ' '*4
				for l in laxis[1::2]: output += '|%-7.1f' % l
			else:
				for l in laxis[0::2]: output += '|%-7.1f' % l
			output += '\n'
		else:
			#head
			output = r'    b\l'
			if laxis[0]%1==0:
				for l in laxis[0::4]: output += ' |%-6.1f' % l
			else:
				output += ' '*2
				for l in laxis[1::4]: output += ' |%-6.1f' % l
			output += '\n'
			#body
			for b,line in zip(baxis,count):
				output += '%+5.1f -' % b if b%1==0 else ' '*7
				for value in line: output += '%2i' % value if value!=0 else ' .'
				output += ' - %+-5.1f\n' % b if b%1!=0 else ' '*8+'\n'
			#rear
			output += ' '*7
			if laxis[0]%1==0:
				for l in laxis[2::4]: output += ' |%-6.1f' % l
			else:
				output += ' '*6
				for l in laxis[3::4]: output += ' |%-6.1f' % l
			output += '\n'

		return output

	def __repr__(self):
		return self.__str__()

	def where(self, all=True):
		cells = []
		if all:
			for idx in np.argwhere(self!=0):
				cells.append('%04i%+04i' % (self.laxis[idx[1]], self.baxis[idx[0]]))
		else:
			lindex, bindex = self.lindex, self.bindex
			laxis = self.laxis[lindex[0]:lindex[1]]
			baxis = self.baxis[bindex[0]:bindex[1]]
			count = self[bindex[0]:bindex[1], lindex[0]:lindex[1]]
			for idx in np.argwhere(count!=0):
				cells.append('%04i%+04i' % (laxis[idx[1]], baxis[idx[0]]))
		return cells


if __name__ == '__main__':
	print('>>check observed cells<<')
	havecube = CellMap('/share/data/mwisp/G???+00/*.fits',l=[0,20])
	print(havecube)

	print('>>change position<<')
	havecube.l = [35,50]
	havecube.b = [-6,6]
	print(havecube)
	
	print('>>search for missing info file in region shown<<')
	haveinfo = CellMap('/share/data/mwisp/infofiles/*info.txt')
	missinginfo = (havecube>0) & (haveinfo==0)
	print(missinginfo)
	print(missinginfo.where(all=False))
