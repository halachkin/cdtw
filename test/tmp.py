from cdtw.src.dtw import*


s = Settings()

s.local_constraint.set_type('p05sym')
s.compute_path = True
s.comute_cost = True
r = [1.23,2.23,3,4]
q = [1.32,2,3.3232,4]


d = dtw(r,q,s)

print d.get_dist()

for i in range(1000):
	d_ = dtw(r,q,s)
	if d.get_dist() != d_.get_dist():
		print "FAIL"
		print d.get_dist()
		print d_.get_dist()


