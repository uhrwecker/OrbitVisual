from env.environment import OrbitalEnvironment
import pygame, sys, numpy

class OrbitVision():
    def __init__(self, display_size=(1000, 800), initial_r=None, initial_v=None, masses=None,
                 seed=123):
        
        numpy.random.seed(seed)
        
        self.width, self.height = display_size
        self.screen = pygame.display.set_mode(display_size)

        r = self.__setup_initial_rs(initial_r)
        v = self.__setup_initial_vs(initial_v)

        self.orbit, self.color, self.masses = self.__setup_orbit(r, v, masses)
        
        self.clock = pygame.time.Clock()
        
    def run(self, nb_stars=20):
        self.star_int = [(numpy.random.randint(0, high=self.width), numpy.random.randint(0, high=self.height))
                         for ind in range(nb_stars)]
        self.game_loop(nb_stars=nb_stars)

    def add_planet(self, r, v, mass):
        self.orbit.add_orbit(r, v, mass)
        self.masses.append(mass)
        self.color = self.color.tolist()
        self.color.append(numpy.random.randint(0, high=255, size=3))        

    def game_loop(self, nb_stars=20):
        
        while True:
    
            for event in pygame.event.get():
                if event.type == pygame.QUIT: sys.exit()

            self.draw_background(number_stars=nb_stars)

            state = self.orbit.step()

            self.draw_planets(state)

            pygame.display.flip()
            self.clock.tick(60)

    def draw_background(self, number_stars=20):
        self.screen.fill((0, 0, 0))
        for index in range(number_stars):
            pygame.draw.circle(self.screen, (255, 255, 255), self.star_int[index], 3)
        
    def draw_planets(self, state):
        for n, st in enumerate(state):
            r, v = st
            radius = 40 + int(str(self.masses[n])[0])*10
            planet = pygame.draw.circle(self.screen, self.color[n], (int(r[0]+self.width/2),
                                                                     int(r[1]+self.height/2)), radius)

    def __setup_initial_rs(self, flag):
        if flag != None:
            return flag
        else:
            r0 = [300, -100, 0]
            r1 = [50, 200, 0]
            r2 = [-300, -300, 0]

            return r0, r1, r2

    def __setup_initial_vs(self, flag):
        if flag != None:
            return flag
        else:
            v0 = [-0.105, 0, 0]
            v1 = [-0.69, 0, 0]
            v2 = [0.5, 0.5, 0]

            return v0, v1, v2

    def __setup_orbit(self, r, v, masses):
        if masses == None:
            masses = [2e27, 1e27, 3e27]
        color = numpy.random.randint(0, high=255, size=(len(masses), 3))
        return OrbitalEnvironment(masses, r, v), color, masses

def main():
    o = OrbitVision()
    o.add_planet([400, 400, 0], [0, 0, 0], 1e10)
    o.run(nb_stars=100)

if __name__ == '__main__':
    main()

