"""
Class interface for stand generation
"""

from math import sqrt

from alinea.adel.stand.stand import agronomicplot

class AgronomicStand(object):
    
    def __init__(self, sowing_density=10, plant_density=10, inter_row=0.8):
        self.sowing_density = sowing_density
        self.inter_row = inter_row
        self.plant_density = plant_density
        self.inter_plant = 1. / inter_row / sowing_density
        
    def plot_dimensions(self, nplants =1, aspect = 'square'):
    
        if aspect =='square':
            side = sqrt(1. / self.sowing_density * nplants)
            nrow = max(1, round(side / self.inter_row))
            plant_per_row = max(1, round(side / self.inter_plant)) 
            plot_length = self.inter_plant * plant_per_row
            plot_width = self.inter_row * nrow
            return plot_length, plot_width
        elif aspect == 'line':
            plot_width = self.inter_row    
            plot_length = (nplants + 1) * self.inter_plant * self.sowing_density / float(self.plant_density) if self.plant_density > 0. else 0.
            return plot_length, plot_width  
        else:
            return 0.5, 0.5
       
    def stand(self, nplants = 1, aspect='square'):
        
        length, width = self.plot_dimensions(nplants, aspect)
        n_emerged, positions, domain, domain_area, convUnit = agronomicplot(length, width, self.sowing_density, self.plant_density, self.inter_row)
        
        return n_emerged, domain, positions, length * width