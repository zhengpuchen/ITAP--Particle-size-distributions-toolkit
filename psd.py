# psd.py
# Class and functions to analyze frequency and cumulative distributions.

import numpy as np
import csv

class Distribution:

    def __init__(self, **kwargs):
        # Get the distribution.
        flag = False
        for key, value in kwargs.items():
            if (key == 'q0'):  # a frequency distribution by number
                self.q = value
                self.by = 'number'
                self.type = 'frequency'
                flag = True
            elif (key == 'q2'):  # a frequency distribution by area
                self.q = value
                self.by = 'area'
                self.type = 'frequency'
                flag = True
            elif (key == 'q3'):  # a frequency distribution by volume
                self.q = value
                self.by = 'volume'
                self.type = 'frequency'
                flag = True
            elif (key == 'Q0'):  # a cumulative distribution by number
                self.Q = value
                self.by = 'number'
                self.type = 'cumulative'
                flag = True
            elif (key == 'Q2'):  # a cumulative distribution by area
                self.Q = value
                self.by = 'area'
                self.type = 'cumulative'
                flag = True
            elif (key == 'Q3'):  # a cumulative distribution by volume
                self.Q = value
                self.by = 'volume'
                self.type = 'cumulative'
                flag = True
        if (flag == False):
            print("Need to specify a distribution type, e.g., q0, q2, q3, Q0, Q2, or Q3")
            exit()

        # Get the sizes corresponding to the distribution.
        flag = False
        for key, value in kwargs.items():
            if (key == 'x'):  # sizes corresponding to the distribution
                self.x = value
                flag = True
        if (flag == False):
            print("Need to specify the sizes corresponding to the distribution, e.g., x.")
            exit()

        # Verify that the size and distribution arrays match.
        if (self.type == 'frequency'):  # This is a frequency distribution.
            if (len(self.x) != len(self.q)+1): 
                print("The frequency distribution and size arrays aren't the correct lengths.")
            else:
                self.x0 = np.empty(len(self.x)-1, dtype=float)
                self.x1 = np.empty(len(self.x0), dtype=float)
                for i in range(len(self.x)-1):
                    self.x0[i] = self.x[i]
                    self.x1[i] = self.x[i+1]
        else:  #  This is a cumulative distribution.
            if (len(self.x) != len(self.Q)):
                print("The cumulative distribution and size arrays aren't the same length.")

        # Verify that the distribution has the correct characteristics.
        if (self.type == 'frequency'):  # We're given a frequency distribution
            integral = self.GetFrequencyDistributionIntegral()
            if (np.abs(integral- 1) > 0.001):
                print("Integral of the frequency distribution isn't equal within 0.1% of one.  Your integral = ", integral)
                exit()
        else:  # We're given a cumulative distribution
            # Verify that the cumulative distribution values are reasonable.
            Q0, Q1, monotonic_flag = self.GetCumulativeDistributionCharacteristics()
            if (Q0 != 0):
                print("The first Q value should equal zero.  Your value is:  ", Q0)
                exit()
            if (Q1 != 1):
                print("The last Q value should equal one.  Your value is:  ", Q1)
                exit()
            if (monotonic_flag == False):
                print("The Q values should increase monotonically.  Yours don't.")
                exit()

        # Generate the opposite type of distribution.
        if (self.type == 'frequency'):  # it's a frequency distribution
            self.GetCumulativeDistributionFromFrequencyDistribution()
        else:  # it's a cumulative distribution
            # Determine the corresponding frequency distribution.
            self.GetFrequencyDistributionFromCumulativeDistribution()
                
        # Determine the average size in the interval.
        self.GetAverageSizeInInterval()

    #####

    def GetAverageSizeInInterval(self, x0=None, x1=None):
        # Determine the average size in the interval.

        if (x0 is not None):  # parameters have been passed
            x_avg = np.empty([len(x0)], dtype=float)
            for i in range(len(x0)):
                x_avg[i] = 0.5*(x0[i] + x1[i])
            return x_avg
        else:
            self.x_avg = np.empty([len(self.x0)], dtype=float)
            for i in range(len(self.x0)):
                self.x_avg[i] = 0.5*(self.x0[i] + self.x1[i])

    #####
    
    def GetFrequencyDistributionIntegral(self, x0=None, x1=None, q=None):
        # Return the integral of the frequency distribution.
        sum = 0
        if (x0 is not None):  # parameters have been passed to the function
            for i in range(len(q)):
                delta_size = x1[i] - x0[i]
                sum = sum + q[i]*delta_size
        else:
            for i in range(len(self.q)):
                delta_size = self.x1[i] - self.x0[i]
                sum = sum + self.q[i]*delta_size
        return sum

    #####

    def GetCumulativeDistributionCharacteristics(self, Q=None):
        # Return the characteristics of the cumulative distribution.
        monotonic_flag = True
        if (Q is not None):  # parameters have been passed to the function
            Q0 = Q[0]
            Q1 = Q[len(Q)-1]
            for i in range(len(Q)-1):
                if (Q[i+1] < Q[i]):
                    monotonic_flag = False 
        else:
            Q0 = self.Q[0]
            Q1 = self.Q[len(self.Q)-1]
            for i in range(len(self.Q)-1):
                if (self.Q[i+1] < self.Q[i]):
                    monotonic_flag = False 
        return Q0, Q1, monotonic_flag

    #####

    def GetD(self, perc):
        perc = perc/100  # convert to a percentage
        
        # Working with the cumulative fraction passing distribution.
        if ((perc < 0) or (perc > 1)):
            raise ValueError ("perc value is outside range [0, 100]")

        i = 0
        while (perc > self.Q[i]):
            i = i+1
        if (i == 0):
            return self.x[0]
        # Linearly interpolate the D value.
        D = self.x[i-1] + (self.x[i]-self.x[i-1])/(self.Q[i]-self.Q[i-1])*(perc-self.Q[i-1])

        return D

    #####
    
    def GetD_Median(self):
        return self.GetD(50)

    #####
    
    def GetSpan(self):
        return ((self.GetD(90) - self.GetD(10))/self.GetD(50))

    #####
    
    def GetD_Mode(self):
        max = 0
        for i in range(len(self.x0)):
            if (self.q[i] > max):
                max = self.q[i]
                imax = i
        return 0.5*(self.x0[imax]+self.x1[imax])

    #####
    
    def GetD_ArithmeticMean(self, x0=None, x1=None, q=None):
        sum = 0
        if (x0 is not None):  # parameters have been passed
            for i in range(len(x0)):
                delta_size = x1[i] - x0[i]
                avg_size = 0.5*(x0[i] + x1[i])
                sum = sum + avg_size*q[i]*delta_size
        else:
            for i in range(len(self.x0)):
                delta_size = self.x1[i] - self.x0[i]
                sum = sum + self.x_avg[i]*self.q[i]*delta_size
        return sum

    #####
        
    def GetD_QuadraticMean(self):
        sum = 0
        for i in range(len(self.x0)):
            delta_size = self.x1[i] - self.x0[i]
            avg_size = 0.5*(self.x0[i] + self.x1[i])
            sum = sum + (avg_size*avg_size)*self.q[i]*delta_size
        return np.sqrt(sum)

    #####
    
    def GetD_CubicMean(self):
        sum = 0
        for i in range(len(self.x0)):
            delta_size = self.x1[i] - self.x0[i]
            avg_size = 0.5*(self.x0[i] + self.x1[i])
            sum = sum + (avg_size*avg_size*avg_size)*self.q[i]*delta_size
        return (sum**(1/3))

    #####
    
    def GetD_HarmonicMean(self):
        sum = 0
        for i in range(len(self.x0)):
            delta_size = self.x1[i] - self.x0[i]
            avg_size = 0.5*(self.x0[i] + self.x1[i])
            sum = sum + (1/avg_size)*self.q[i]*delta_size
        return (1/sum)

    #####
        
    def GetD_GeometricMean(self):
        sum = 0
        for i in range(len(self.x0)):
            delta_size = self.x1[i] - self.x0[i]
            avg_size = 0.5*(self.x0[i] + self.x1[i])
            sum = sum + np.log(avg_size)*self.q[i]*delta_size
        return (np.exp(sum))

    #####
        
    def GetStandardDeviation(self):
        # Determine the standard deviation from the frequency distribution.
        
        xbar = self.GetD_ArithmeticMean()
        sum = 0
        for i in range(len(self.x0)):
            delta_size = self.x1[i] - self.x0[i]  # size of this bin
            avg_size = self.x_avg[i]
            sum = sum + ((avg_size-xbar)**2)*self.q[i]*delta_size
        return np.sqrt(sum)

    #####
    
    def GetSkewness(self):

        xbar = self.GetD_ArithmeticMean()
        stdev = self.GetStandardDeviation()
        sum = 0
        for i in range(len(self.x0)):
            delta_size = self.x1[i] - self.x0[i]
            avg_size = self.x_avg[i]
            sum = sum + ((avg_size-xbar)**3)*self.q[i]*delta_size
        return (sum/(stdev**3))

    #####
        
    def GetExcessKurtosis(self):

        xbar = self.GetD_ArithmeticMean()
        stdev = self.GetStandardDeviation()
        sum = 0
        for i in range(len(self.x0)):
            delta_size = self.x1[i] - self.x0[i]
            avg_size = self.x_avg[i]
            sum = sum + ((avg_size-xbar)**4)*self.q[i]*delta_size
        kurtosis = sum/(stdev**4)
        return (kurtosis-3)

    #####

    def GetFrequencyDistributionFromCumulativeDistribution(self, x=None, Q=None):

        if (x is not None):  # parameters have been passed to the function
            # First check that the cumulative distribution characteristics are reasonable.
            Q0, Q1, monotonic_flag = self.GetCumulativeDistributionCharacteristics(Q)
            # Verify that the cumulative distribution values are reasonable.
            if (Q0 != 0):
                raise ValueError ("The first Q value should equal zero.  Your value is:  ", Q0)
            if (Q1 != 1):
                raise ValueError ("The last Q value should equal one.  Your value is:  ", Q1)
            if (monotonic_flag == False):
                raise ValueError ("The Q values should increase monotonically.  Yours don't.")

            # Determine the corresponding frequency distribution.
            x0 = np.empty([len(x)-1], dtype=float)
            x1 = np.empty([len(x0)], dtype=float)
            q = np.empty([len(x0)], dtype=float)

            for i in range(len(x)-1):
                x0[i] = x[i]
                x1[i] = x[i+1]
                delta_size = x1[i] - x0[i]
                q[i] = (Q[i+1] - Q[i])/delta_size
            return x0, x1, q
        else:
            self.x0 = np.empty([len(self.x)-1], dtype=float)
            self.x1 = np.empty([len(self.x0)], dtype=float)
            self.q = np.empty([len(self.x0)], dtype=float)
            
            # First check that the cumulative distribution characteristics are reasonable.
            Q0, Q1, monotonic_flag = self.GetCumulativeDistributionCharacteristics()
            if (Q0 != 0):
                raise ValueError ("The first Q value should equal zero.  Your value is:  ", Q0)
            if (Q1 != 1):
                raise ValueError ("The last Q value should equal one.  Your value is:  ", Q1)
            if (monotonic_flag == False):
                raise ValueError ("The Q values should increase monotonically.  Yours don't.")
            
            for i in range(len(self.x)-1):
                self.x0[i] = self.x[i]
                self.x1[i] = self.x[i+1]
                delta_size = self.x1[i] - self.x0[i]
                self.q[i] = (self.Q[i+1] - self.Q[i])/delta_size

    #####

    def GetCumulativeDistributionFromFrequencyDistribution(self, x=None, q=None):

        if (x is not None):  # parameters have been passed to the
            # function First check that the frequency distribution
            # characteristics are reasonable.

            integral = self.GetFrequencyDistributionIntegral()
            if (np.abs(integral- 1) > 0.001):
                print("Integral of the frequency distribution isn't equal within 0.1% of one.  Your integral = ", integral)
                exit()
            
            # Determine the corresponding cumulative distribution.
            Q = np.empty([len(x)], dtype=float)

            for i in range(len(x0)):
                x[i] = x0[i]
                if (i > 0):
                    delta_size = x1[i-1] - x0[i-1]
                    Q[i] = Q[i-1] + q[i-1]*delta_size
                else:
                    Q[i] = 0
            i = len(x)-1
            x[i] = x1[i-1]
            delta_size = x1[i-1] - x0[i-1]
            Q[i] = Q[i-1] + q[i-1]*delta_size
            return x, Q
        else:
            self.Q = np.empty([len(self.x)], dtype=float)

            for i in range(len(self.x0)):
                self.x[i] = self.x0[i]
                if (i > 0):
                    delta_size = self.x1[i-1] - self.x0[i-1]
                    self.Q[i] = self.Q[i-1] + self.q[i-1]*delta_size
                else:
                    self.Q[i] = 0
            i = len(self.x)-1
            self.x[i] = self.x1[i-1]
            delta_size = self.x1[i-1] - self.x0[i-1]
            self.Q[i] = self.Q[i-1] + self.q[i-1]*delta_size
                
    #####

    def GetFrequencyDistributionByNumber(self):
        if (self.by == 'number'):
            print("Already have a number distribution.")
            return self.q
        elif (self.by == 'area'):
            print("Converting an area distribution to a number distribution.")
            print("Assuming the particle shape doesn't change across size.")
            # Convert an area distribution to a number distribution.
            q0 = np.empty(len(self.q), dtype=float)

            # First calculate the denominator.
            sum = 0
            for i in range(len(self.q)):
                deltax = self.x1[i] - self.x0[i]
                sum = sum + self.q[i]*deltax/(self.x_avg[i]**2)

            # Now calculate the distribution.
            for i in range(len(self.q)):
                q0[i] = (self.q[i]/(self.x_avg[i]**2))/sum
            return q0
        else:
            # Convert a volume distribution to a number distribution.
            print("Converting a volume distribution to a number distribution.")
            print("Assuming the particle shape doesn't change across size.")
            q0 = np.empty(len(self.q), dtype=float)

            # First calculate the denominator.
            sum = 0
            for i in range(len(self.q)):
                deltax = self.x1[i] - self.x0[i]
                sum = sum + self.q[i]*deltax/(self.x_avg[i]**3)

            # Now calculate the distribution.
            for i in range(len(self.q)):
                q0[i] = self.q[i]/(self.x_avg[i]**3)/sum

            return q0

    #####
    
    def GetFrequencyDistributionByArea(self):
        if (self.by == 'area'):
            print("Already have an area distribution.")
            return self.q
        elif (self.by == 'number'):
            print("Converting a number distribution to an area distribution.")
            print("Assuming the particle shape doesn't change across size.")
            # Convert a number distribution to an area distribution.
            q2 = np.empty(len(self.q), dtype=float)

            # First calculate the denominator.
            sum = 0
            for i in range(len(self.q)):
                deltax = self.x1[i] - self.x0[i]
                sum = sum + self.q[i]*deltax*(self.x_avg[i]**2)

            # Now calculate the distribution.
            for i in range(len(self.q)):
                q2[i] = (self.q[i]*(self.x_avg[i]**2))/sum
            return q2
        else:
            # Convert a volume distribution to an area distribution.
            print("Converting a volume distribution to an area distribution.")
            print("Assuming the particle shape doesn't change across size.")
            q2 = np.empty(len(self.q), dtype=float)

            # First calculate the denominator.
            sum = 0
            for i in range(len(self.q)):
                deltax = self.x1[i] - self.x0[i]
                sum = sum + self.q[i]*deltax/self.x_avg[i]

            # Now calculate the distribution.
            for i in range(len(self.q)):
                q2[i] = (self.q[i]/self.x_avg[i])/sum
            return q2

    #####
            
    def GetFrequencyDistributionByVolume(self):
        if (self.by == 'volume'):
            print("Already have a volume distribution.")
            return self.q
        elif (self.by == 'number'):
            # Convert a number distribution to a volume distribution.
            print("Assuming the particle shape doesn't change across size.")
            q3 = np.empty(len(self.q), dtype=float)

            # First calculate the denominator.
            sum = 0
            for i in range(len(self.q)):
                deltax = self.x1[i] - self.x0[i]
                sum = sum + (self.x_avg[i]**3)*self.q[i]*deltax

            # Now calculate the distribution.
            for i in range(len(self.q)):
                q3[i] = ((self.x_avg[i]**3)*self.q[i])/sum
            return q3
        else:
            # Convert an area distribution to a volume distribution.
            print("Assuming the particle shape doesn't change across size.")
            q3 = np.empty(len(self.q), dtype=float)
        
            # First calculate the denominator.
            sum = 0
            for i in range(len(self.q)):
                deltax = self.x1[i] - self.x0[i]
                sum = sum + (self.x_avg[i])*self.q[i]*deltax

            # Now calculate the distribution.
            for i in range(len(self.q)):
                q3[i] = self.x_avg[i]*self.q[i]/sum
            return q3

    #####

    def WriteFrequencyDistribution(self, filename=None):
        if (filename is None):  # no parameter has been passed
            filename = "FrequencyDistribution.csv"
        with open(filename, 'w', newline='') as csv_file:
            print("Writing frequency distribution to:  ", filename)
            csv_writer = csv.writer(csv_file, delimiter=',')
            for i in range(len(self.q)):
                csv_writer.writerow([self.x0[i], self.x1[i], self.q[i]])

    #####
                
    def WriteCumulativeDistribution(self, filename=None):
        if (filename is None):  # no parameter has been passed
            filename = "CumulativeDistribution.csv"
        with open(filename, 'w', newline='') as csv_file:
            print("Writing cumulative distribution to:  ", filename)
            csv_writer = csv.writer(csv_file, delimiter=',')
            for i in range(len(self.Q)):
                csv_writer.writerow([self.x[i], self.Q[i]])

    #####

    def GetSauterMeanDiameter(self):
        q2 = self.GetFrequencyDistributionByArea()
        return self.GetD_ArithmeticMean(self.x0, self.x1, q2)

    #####

    def GetDeBrouckereMeanDiameter(self):
        q3 = self.GetFrequencyDistributionByVolume()
        return self.GetD_ArithmeticMean(self.x0, self.x1, q3)

    #####

    def GetDistributionMoment(self, k, r=None):

        if (r is not None):
            if (r not in [0, 2, 3]):
                raise ValueError ("r should be 0, 2, or 3.  You have: ", r)
            elif (r == 0):
                qr = self.GetFrequencyDistributionByNumber()
            elif (r == 2):
                qr = self.GetFrequencyDistributionByArea()
            else:
                qr = self.GetFrequencyDistributionByVolume()
            
            sum = 0
            for i in range(len(qr)):
                deltax = self.x1[i] - self.x0[i]
                sum = sum + (self.x_avg[i]**k)*qr[i]*deltax
            return sum
        else:
            sum = 0
            for i in range(len(self.q)):
                deltax = self.x1[i] - self.x0[i]
                sum = sum + (self.x_avg[i]**k)*self.q[i]*deltax
            return sum

    #####

    def FractionInInterval(self, x_lo, x_hi):
        # Find the fraction of values in the range [x_lo, x_hi].

        if (x_hi < x_lo):  
            print("Error in:  GetFractionInInterval:  x_hi < x_lo")
            exit(1)
            
        i = 0
        while ((x_lo > self.x[i]) and (i < len(self.x)-1)):
            i = i + 1
        if (i == 0):
            Q_lo = self.Q[0]
        elif (i == len(self.x)-1):
            #  The x_lo value is beyond the end of the distribution.
            return 0
        else:
            # Linearly interpolate the Q_lo value.
            Q_lo = self.Q[i-1] + (self.Q[i]-self.Q[i-1])/(self.x[i]-self.x[i-1])*(x_lo-self.x[i-1])

        i = 0
        while ((x_hi > self.x[i]) and (i < len(self.x)-1)):
            i = i + 1
        if (i == 0):
            Q_hi = self.Q[0]
        elif (i == len(self.x)-1):
            #  The x_hi value is beyond the end of the distribution.
            Q_hi = self.Q[len(self.x)-1]
        else:
            # Linearly interpolate the Q_hi value.
            Q_hi = self.Q[i-1] + (self.Q[i]-self.Q[i-1])/(self.x[i]-self.x[i-1])*(x_hi-self.x[i-1])
        return (Q_hi - Q_lo)
