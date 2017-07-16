#reference values are from DTW - R projekt :
#http://dtw.r-forge.r-project.org/


import numpy as np

from cdtw.cydtw import Settings
from cdtw.pydtw import *

class TestSuite:
    """
    Create a suite of tests 
    """
    def __init__(self):
        """
        Creates a test suite object
        """
        self.tests = 0
        self.failures = 0   
    
    def run_test(self, computed, expected, message = ""):
        """
        compare computed and expected values
        """
        def round_inputs(comp,exp):
            if isinstance(exp, (int, long, float, complex)):
                return round(comp,2),round(exp,2)
            elif type(exp) == list and type(exp[0]) != list:
                return comp,exp
            elif type(exp[0]) == list:
                return [[round(elem,2) for elem in row] for row in comp], [[round(elem,2) for elem in row] for row in exp]
            else:
                return comp,exp

        computed_, expected_ = round_inputs(computed,expected)
        self.tests += 1
        if computed_ != expected_:
            print message + " Computed: " + str(computed) + \
                            " Expected: " + str(expected)
            self.failures += 1
    
    def report(self):
        """
        Report results
        """
        print "Ran " + str(self.tests) + " tests. " \
                     + str(self.failures) + " failures."
    def get_failures_count(self):
        """
        failures count getter
        """
        return self.failures


def run(settings, test_suite):
    ref =   [2.323, 3.555, 1.454, 2.579, 1.214, 6.214] 
    query = [1.205, 1.921, 7.454, 3.214, 4.104, 2.222]

    #true dist, local constraint - symmetric1 (dp2)
    dist = 12.434
    pathx = np.array([1, 1, 2, 2, 2, 3, 4, 5, 6]) - 1
    pathy = np.array([1, 2, 3, 4, 5, 6, 6, 6, 6]) - 1
    path = zip(pathx,pathy)

    s.step.set_type('dp2')

    

    d = dtw(ref,query,s)

    dist_ = d.get_dist()
    path_ = d.get_path()

    test_id = 'symmetric1 (dp2)'
    t.run_test(round(dist_,3), dist, str("Distance error, test ID: " + str(test_id) ))
    t.run_test(path_, path, str("Path error, test ID: " + str(test_id)))


    #true dist, local constraint - Sakoe Chiba P = 0 sym
    dist = 14.775
    # here is the path from R package
    pathx = [0, 0, 0, 1, 1, 2, 3, 4, 5]
    pathy = [0, 1, 2, 3, 4, 5, 5, 5, 5]
    path = zip(pathx,pathy)


    #computing values
    s.step.set_type('p0sym')

    d = dtw(ref,query,s)
    dist_ = d.get_dist()
    path_ = d.get_path()
    np.set_printoptions(linewidth = 180, suppress = True, precision  = 3)
    

    #tests
    test_id = 'Sakoe-Chiba P = 0 sym'
    t.run_test(round(dist_,3), dist, str("Distance error, test ID: " + str(test_id) ))
    t.run_test(path_, path, str("Path error, test ID: " + str(test_id)))



    #local constraint - Sakoe Chiba P = 0 asym
    #true dist
    dist = 5.824
    # here is the path from R package, which is probably wrong
    pathx = np.array([1,2,3,4,5,6,6,6,6]) - 1
    pathy = np.array([1,2,2,2,2,3,4,5,6]) - 1
    path = zip(pathx,pathy)

    #computing values
    s.step.set_type('p0asym')
    
    d = dtw(ref,query,s)
    dist_ = d.get_dist()
    path_ = d.get_path()

    # 

    #tests
    test_id = 'Sakoe-Chiba P = 0 asym'
    t.run_test(round(dist_,3), dist, str("Distance error, test ID: " + str(test_id) ))
    t.run_test(path_, path, str("Path error, test ID: " + str(test_id)))


        #true dist#local constraint - Sakoe Chiba P = 1/2 sym
    dist = 19.640
    pathx = np.array([1,2,2,2,3,4,5,6]) - 1 
    pathy = np.array([1,2,3,4,5,6,6,6]) - 1
    path = zip(pathx,pathy)

    #computing values
    s.step.set_type('p05sym')
    
    d = dtw(ref,query,s)
    dist_ = d.get_dist()
    path_ = d.get_path()


    #tests
    test_id = 'Sakoe-Chiba P = 1/2 sym'
    t.run_test(round(dist_,3), dist, str("Distance error, test ID: " + str(test_id) ))
    t.run_test(path_, path, str("Path error, test ID: " + str(test_id)))


    #local constraint - Sakoe Chiba P = 1/2 asym
    #true dist
    dist = 10.564
    pathx = np.array([1,2,3,4,4,4,5,6]) - 1
    pathy = np.array([1,2,2,3,4,5,6,6]) - 1

    path = zip(pathx,pathy)

    #computing values
    s.step.set_type('p05asym')
    
    d = dtw(ref,query,s)
    dist_ = d.get_dist()
    path_ = d.get_path()
    
    #tests
    test_id = 'Sakoe-Chiba P = 1/2 asym'
    t.run_test(round(dist_,3), dist, str("Distance error, test ID: " + str(test_id) ))
    t.run_test(path_, path, str("Path error, test ID: " + str(test_id)))


    #local constraint - Sakoe Chiba P = 1 sym
    #true dist
    dist = 20.863

    pathx = np.array([1,2,2,3,4,5,6]) - 1
    pathy = np.array([1,2,3,4,5,6,6]) - 1
    path = zip(pathx,pathy)

    #computing values
    s.step.set_type('p1sym')
    
    d = dtw(ref,query,s)
    dist_ = d.get_dist()
    path_ = d.get_path()
    #tests
    test_id = 'Sakoe-Chiba P = 1 sym'
    t.run_test(round(dist_,3), dist, str("Distance error, test ID: " + str(test_id) ))
    t.run_test(path_, path, str("Path error, test ID: " + str(test_id)))



    #local constraint - Sakoe Chiba P = 1 asym

    pathx = np.array([1,2,2,3,4,5,6]) - 1
    pathy = np.array([1,2,3,4,5,6,6]) - 1
    #true dist
    dist = 12.1695
    path = zip(pathx,pathy)

    #computing values
    s.step.set_type('p1asym')
    
    d = dtw(ref,query,s)
    dist_ = d.get_dist()
    path_ = d.get_path()

    #tests
    test_id = 'Sakoe-Chiba P = 1 asym'
    t.run_test(round(dist_,4), dist, str("Distance error, test ID: " + str(test_id) ))
    t.run_test(path_, path, str("Path error, test ID: " + str(test_id)))



    #local constraint - Sakoe Chiba P = 2 sym
    #look like R package is wrong in this case

    #true dist
    dist = 27.204
    pathx = np.array([1,2,3,3,4,5,6]) - 1 
    pathy = np.array([1,2,3,4,5,6,6]) - 1
    path = zip(pathx,pathy)

    #computing values
    s.step.set_type('p2sym')
    
    d = dtw(ref,query,s)
    dist_ = d.get_dist()
    path_ = d.get_path()

    #tests
    test_id = 'Sakoe-Chiba P = 2 sym'
    t.run_test(round(dist_,3), dist, str("Distance error, test ID: " + str(test_id) ))
    t.run_test(path_, path, str("Path error, test ID: " + str(test_id)))




    #local constraint - Sakoe Chiba P = 2 asym
    #true dist
    dist = 13.90567
    pathx = np.array([1,2,3,3,4,5,6]) - 1
    pathy = np.array([1,2,3,4,5,6,6]) - 1
    path = zip(pathx,pathy)

    #computing values
    s.step.set_type('p2asym')
    
    d = dtw(ref,query,s)
    dist_ = d.get_dist()
    path_ = d.get_path()

    # print np.array(d.get_cost())
    #tests
    test_id = 'Sakoe-Chiba P = 2 asym'
    t.run_test(round(dist_,5), dist, str("Distance error, test ID: " + str(test_id) ))
    t.run_test(path_, path, str("Path error, test ID: " + str(test_id)))

    np.set_printoptions(linewidth = 180, suppress = True, precision  = 3)

    ref = [1,0,700,800,600,451,0,0,0,0,100]
    query = [0,0,0,0,0,150,200,250,100,10,10]
    dist = 4082
    pathx = np.array([1, 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9, 10, 10, 10, 11, 11, 11]) - 1
    pathy = np.array([1, 2 , 2 , 2 , 3 , 3 , 3 , 4 , 5,  6,  7,  8,  9, 10, 11]) - 1
    path = zip(pathx,pathy)

    s.step.set_type('p05sym')
    d = dtw(ref,query,s)
    dist_ = d.get_dist()
    path_ = d.get_path()

    # print np.array(d.get_cost())
    #tests
    test_id = 'Sakoe-Chiba P = 0.5 sym'
    t.run_test(round(dist_,5), dist, str("Distance error, test ID: " + str(test_id) ))
    t.run_test(path_, path, str("Path error, test ID: " + str(test_id)))

    t.report()
   

s = Settings()
t = TestSuite()

s.compute_path = True
run(s,t)



 