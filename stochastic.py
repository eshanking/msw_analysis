#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 12:28:02 2022

@author: beckettpierce
"""
import numpy as np
import matplotlib.pyplot as plt
import random
class population():
    def __init__(self,
                     pop_count = None,
                     a_count = None,
                     b_count = None,
                     total_time = None,
                     time_step = None,
                     rel_fitness = None):
        self.pop_count = pop_count
        self.a_count = a_count
        self.total_time = total_time
        self.time_step = time_step
        self.rel_fitness = rel_fitness
        self.b_count = b_count
    def one_step(self):
        #need to calculate birth and death of a & b chance
        prob_birth_a = ((self.a_count*self.rel_fitness)*(self.pop_count-self.a_count))/((self.rel_fitness*self.a_count+self.pop_count-self.a_count)*(self.pop_count))
        prob_birth_b = 1-prob_birth_a
        prob_death_a = ((self.pop_count-self.a_count)*(self.a_count)) / ((self.rel_fitness*self.a_count+self.pop_count-self.a_count)*(self.pop_count))
        prob_death_b = 1-prob_death_a
        chooser = np.random.random(1)[0]
        chooser2 = np.random.random(1)[0]
        if chooser <= prob_birth_a:
            self.a_count = self.a_count+1
            self.b_count = self.b_count-1
        if chooser <= prob_death_a:
            self.a_count = self.a_count-1
            self.b_count = self.b_count+1
        return self.a_count,self.b_count
    def run_sim(self):
        time=0
        sim_results = np.zeros((int(np.floor(self.total_time/self.time_step)),3))
        counter = 0
        while time < self.total_time:
            a_count, b_count = self.one_step()
            sim_results[counter,0] = time
            sim_results[counter,1] = a_count
            sim_results[counter,2] = b_count
            
            time = time+ self.time_step
            counter = counter + 1
        return sim_results
example_sim = population(pop_count=1000,
                         a_count=1,
                         b_count=999,
                         rel_fitness=1.5,
                         time_step=1,
                         total_time=100000
                         )
results  = example_sim.run_sim()