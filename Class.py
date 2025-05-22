"""
Tumor-Immune Interaction Simulation Model

This code simulates the interaction between tumor cells and the immune system,
focusing on cancer cell growth, mutation, and the immune response through
B cells and antibodies in a 3D spatial environment.

Code is by no means the most efficient, but it works.
"""

import math
import numpy as np
import random
import itertools


# Cell Classes

class Cancer:
    """
    Represents a cancer cell with specific antigenic features and replication properties.
    
    Cancer cells can replicate, mutate, and potentially die due to the presence of antibodies.
    """
    species = "c"  # Identifier for cancer cells
    
    def __init__(self, Feature, mut_p, rep):
        """
        Initialize a cancer cell with specific features and properties.
        
        Parameters:
        - Feature: Feature object representing the cell's antigenic profile
        - mut_p: Mutation probability during replication
        - rep: Base replication rate
        """
        self.Feature = Feature  # Antigenic features
        self.replication_rate = rep
        self.mut_p = mut_p  # Mutation probability
    
    def replicate(self):
        """
        Create a new cancer cell through replication, potentially with mutations.
        
        The mutation process can alter the antigenic features of the daughter cell,
        potentially allowing it to escape immune detection or enter different zones.
        
        Returns:
        - A new Cancer cell object, potentially with mutated features
        """
        # Scale mutation probability by square root of dimensions to maintain consistent average mutation distance regardless of feature dimension
        mutp = self.mut_p / math.sqrt(self.Feature.d)
        
        # Generate mutation mask based on probability
        mutate_mask = np.random.rand(len(self.Feature.features))
        
        # Apply mutations to features based on mask
        # If mutate < mutp, apply a mutation scaled by mutp; otherwise keep the original value
        new_f = [i + (mutate/mutp - 0.5) if mutate < mutp else i for i, mutate in zip(self.Feature.features, mutate_mask)]
        
        # Return new cancer cell with the features
        return Cancer(Feature(new_f), self.mut_p, self.replication_rate)
    
    def die(self, grid, Grids):
        """
        Handle the death of a cancer cell.
        
        Parameters:
        - grid: The Grid object containing this cancer cell
        - Grids: The Grids object managing the entire physical space
        """
        # Remove this grid from the list of cancer locations
        Grids.cancer_list.remove(grid)
        # Decrement cancer cell count
        Grids.cancer_n -= 1
        # Remove the cell from the grid
        grid.cell = None

class Bcell:
    """
    Represents a B cell in the immune system that can recognize cancer antigens.
    
    In fact Bcell plays no part at all by itself and we only needed the Antibody class. But we had it originally so it is still here. 
    """
    species = "b"  # Identifier for B cells
    
    def __init__(self, Feature):
        """
        Initialize a B cell with specific features and properties.
        
        Parameters:
        - Feature: Feature object representing the cell's recognition profile
        """
        self.Feature = Feature  # Antigenic recognition profile
    

class Antibody:
    """
    Represents an antibody produced by B cells that can bind to and mark cancer cells.
    
    Antibodies have a specific recognition profile, effectiveness, and lifecycle states
    (inactive and active).
    """
    def __init__(self, Feature, time, state, effectiveness, multiplier=1):
        """
        Initialize an antibody with specific features and properties.
        
        Parameters:
        - Feature: Feature object representing the antibody's binding profile
        - time: Initial lifetime in current state
        - state: Initial state ("active" or "inactive")
        - effectiveness: Binding effectiveness parameter
        - multiplier: Effectiveness multiplier (default 1)
        """
        self.state = state  # "active" or "inactive"
        self.Feature = Feature  # Binding profile
        self.life_time = time  # Time remaining in current state
        self.effectiveness = effectiveness  # Binding effectiveness
        self.multiplier = multiplier  # Effectiveness multiplier
    
    def activate(self):
        """
        Activate the antibody, changing its state and setting its lifetime.
        """
        self.state = "active"
    
    
    def flow(self):
        """
        Handle the antibody's flow in the bloodstream while inactive.
        
        If lifetime (remaining time of inactivity) reaches 0, the antibody is activated.
        """
        if self.life_time == 0:
            self.activate()
        else:
            self.life_time -= 1
    
    def update(self):
        """
        Update the antibody's state based on its current state.
        """
        if self.state == "inactive":
            # Handle inactive state - continue flowing in bloodstream
            self.flow()
    
    
# Feature Classes
        
class Feature:
    """
    Represents a point in the feature space, defining antigenic properties of cells or antibodies.
    
    Features can be thought of as the molecular characteristics that determine
    how cells interact with each other and with the immune system.
    """
    def __init__(self, features):
        """
        Initialize a Feature object with specific coordinates in feature space.
        
        Parameters:
        - features: Tuple or list of feature values representing coordinates in feature space
        """
        self.features = features
        self.d = len(features)  # Dimension of the space
    
    def distance(self, other):
        """
        Calculate the L1 (Manhattan) distance between this Feature and another Feature.
        
        This distance represents how different far two antigenic profiles are from each other,
        which affects immune recognition.
        
        Parameters:
        - other: Another Feature object to compare with
        
        Returns:
        - The L1 distance between the two Features
        """
        return np.linalg.norm(np.array(self.features) - np.array(other.features), 1)


# Grid System Classes and Functions
# This is about the physical site where cancer cells reside

def locate_space(Grids, cord):
    """
    Find an empty neighboring grid cell for a new cell to occupy.
    
    This function is used when cells replicate to find space for the daughter cell.
    
    Parameters:
    - Grids: The Grids object managing the entire simulation space
    - cord: Coordinates (x,y,z) of the parent cell
    
    Returns:
    - Coordinates of an empty neighboring grid cell, or None if the randomly chosen neighbor isn't empty or valid
    """
    x, y, z = cord
    # Choose a random direction in 3D space
    x_, y_, z_ = (random.randint(-1, 1), random.randint(-1, 1), random.randint(-1, 1))
    nx, ny, nz = (x + x_, y + y_, z + z_)
    
    # Check if the new coordinates are valid and the cell is empty
    if Grids.valid_cord((nx, ny, nz)): # Valid as in, not outside the entire space, which is a finite n by n by n space.
        if Grids.grids[nx][ny][nz].cell is None:
            return nx, ny, nz  # Return coordinates of empty neighboring cell

class Grid:
    """
    Represents a single unit volume in the 3D physical space.
    
    Each Grid can contain a cancer cell and up to roughly 1000 antibodies, and tracks its own time.
    """
    def __init__(self, cord, antibody_dens=[], cell=None):
        """
        Initialize a Grid at specific coordinates.
        
        Parameters:
        - cord: Coordinates (x,y,z) of this grid in the simulation space
        - antibody_dens: Initial list of antibodies in this grid (default empty), not used
        - cell: Initial cell occupying this grid (default None)
        """
        self.cord = cord  # Coordinates
        self.time = 0  # Local time counter
        self.antibody = [x for x in antibody_dens]  # List of antibodies
        self.cell = cell  # Cell occupying this grid (if any)
    
    def update_time(self):
        """
        Increment the local time counter. To make sure a newly born cancer cell don't replicate again in the same turn
        """
        self.time += 1
    
    def add_antibody(self, antibodies):
        """
        Add antibodies to this grid.
        
        Parameters:
        - antibodies: List of Antibody objects to add
        """
        if len(self.antibody)>1e3:
            self.antibody=self.antibody[100:]  # Prevent  memory overflow, and too much antibodies cooexisting.
        self.antibody += antibodies
    
    def antibody_effectiveness(self):
        """
        Calculate the effectiveness of antibodies against the cell in this grid. This effectiveness is the death probability of the cancer cell
        
        Returns:
        - Tuple of (effectiveness probability, index of last antibody considered)
        - (-1, 0) if no cell or no antibodies are present
        """
        # If no cell or no antibodies, return default values
        if not(self.cell and self.antibody):
            return -1, 0
        
        effectiveness = 0
        for body, i in zip(self.antibody, range(len(self.antibody))):
            if body.state == "active":
                # Calculate distance between antibody and cell features
                distance = self.cell.Feature.distance(body.Feature)
                # Calculate effectiveness based on distance and antibody properties
                effectiveness += body.multiplier * math.exp(-body.effectiveness * distance)
                
                # If effectiveness/death prob reaches 1, return 1
                if effectiveness >= 1:
                    return 1, i
        
        # Return calculated effectiveness and index of last antibody
        return effectiveness, i
    
    def cancer_action(self, Grids):
        """
        Handle actions related to cancer cells in this grid, including potential death by antibodies.
        
        Parameters:
        - Grids: The Grids object managing the entire simulation space
        """
        # Calculate antibody effectiveness against the cancer cell
        prob, used = self.antibody_effectiveness()
        grids = Grids.grids
        x, y, z = self.cord
        grid = grids[x][y][z]
        
        # Randomly determine if cell dies based on antibody effectiveness
        if random.random() < prob:
            self.cell.die(self, Grids)
            # Remove used antibodies
            self.antibody = self.antibody[used:]
            return
    
    def replicate_cell(self, cell, Grids):
        """
        Handle cell replication, creating a daughter cell in a neighboring grid.
        
        Parameters:
        - cell: The cell that is replicating
        - Grids: The Grids object managing the entire simulation space
        """
        # Create new cell through replication (potentially with mutations)
        new_cell = self.cell.replicate()
        
        # Find space for the new cell
        space = locate_space(Grids, self.cord)
        if space:
            x, y, z = space
            # Place new cell in the found space
            Grids.grids[x][y][z].cell = new_cell
            # Update time to prevent double update in the same round
            Grids.grids[x][y][z].time += 1
            
            # Update counters based on cell type
            if cell.species == "c":  # Cancer cell
                Grids.cancer_n += 1
                Grids.cancer_list.append(Grids.grids[x][y][z])
    
    def update(self, Grids):
        """
        Update this grid for the current time step, handling cell actions.
        
        Parameters:
        - Grids: The Grids object managing the entire simulation space
        """
        time = Grids.time
        # Skip update if the cancer cell here is a new born
        if self.time > time:
            return
        
        # Handle cancer cell actions if a cell is present
        if self.cell:
            # Check for cancer cell replication based on replication rate
            if random.random() < self.cell.replication_rate:
                self.replicate_cell(self.cell, Grids)
            
            # Handle cancer-specific actions
            if self.cell.species == "c":
                self.cancer_action(Grids)
        
        # Update local time
        self.update_time()
    
    def update_anti(self):
        """
        Update antibodies in this grid, handling their lifecycle and cleanup.
        
        Parameters:
        - Grids: The Grids object managing the entire simulation space
        """
        if self.antibody:
            # Update each antibody, decreasing their self inactive timer if they are inactive
            for body in self.antibody:
                body.update()

class Grids:
    """
    Manages the entire 3D simulation space, containing all Grid objects and global state.
    
    This class handles the overall simulation, tracking cancer cells, antibodies,
    and other global statistics.
    """
    def __init__(self, n, cancer_cord, cancer_feature, mut_p, c_rep):
        """
        Initialize the simulation space with initial cancer cells.
        
        Parameters:
        - n: Size of the cubic simulation space (nxnxn)
        - cancer_cord: List of coordinates for initial batch of cancer cells
        - cancer_feature: homogenous Feature object for initial cancer cells, same for all of the initial batch
        - mut_p: Mutation probability for cancer cells
        - c_rep: replication rate for cancer cells
        """
        # Create initial cancer cell template
        initial_cancer = Cancer(cancer_feature, mut_p, c_rep)
        
        # Initialize simulation parameters
        self.n = n  # Grid size
        self.time = 0  # Global time counter

        
        # Create 3D grid of Grid objects
        self.grids = np.empty((n, n, n), dtype=Grid)
        # Create a virtual grid to keep track of all the antibodies, where they aren't removed with cancer cells
        self.virtual_grid = Grid((-1, -1, -1))
        
        # Initialize all grid cells
        for x in range(n):
            for y in range(n):
                for z in range(n):
                    self.grids[x][y][z] = Grid((x, y, z))
        
        # Initialize cancer tracking
        self.cancer_n = len(cancer_cord)  # Number of cancer cells
        self.cancer_list = []  # List of grids containing cancer cells
        
        # Place initial cancer cells
        for cord in cancer_cord:
            x, y, z = cord
            self.grids[x][y][z].cell = initial_cancer
            self.cancer_list.append(self.grids[x][y][z])
        
        # Create flat list of all grids for easy iteration
        self.grid_list = []
        for x in range(self.n):
            for y in range(self.n):
                for z in range(self.n):
                    self.grid_list.append(self.grids[x][y][z])
        
        
        self.antigen_sampling=False #To start, the grid is not yet sampling the antigens
    def sample_antibody_dist(self):
        """
        Sample the antibody distribution across the simulation space. This is just for visualization and does not affect simulation
        
        Returns:
        - List of antibodies from sampled grid cells
        """
        # Sample indices for representative grid cells
        indices = [1, 500, 1000, 1500, -1500, -1000, -500, -1] 
        sample = []
        
        # Collect antibodies from sampled grids
        for i in indices:
            grid = self.grid_list[i]
            sample += grid.antibody
        
        return sample
    
    def cancer_distribution(self):
        """
        Get the feature distribution of all cancer cells.
        
        Returns:
        - List of feature vectors for all cancer cells
        """
        features = []
        for grid in self.get_cancer_list():
            features.append(grid.cell.Feature.features)
        return features
    
    def get_list(self):
        """
        Get the flat list of all grid cells.
        
        Returns:
        - List of all Grid objects
        """
        return self.grid_list
    
    def get_cancer_list(self):
        """
        Get the list of grid cells containing cancer cells.
        
        Returns:
        - List of Grid objects containing cancer cells
        """
        return self.cancer_list
    
    def valid_cord(self, cord):
        """
        Check if coordinates are valid within the simulation space.
        
        Parameters:
        - cord: Coordinates (x,y,z) to check
        
        Returns:
        - Boolean indicating whether coordinates are valid
        """
        for c in cord:
            if not (0 <= c < self.n):
                return False
        return True
    
    def sample_antigens(self):
        """
        Sample antigens from cancer cells for immune recognition.
        
        This simulates the process of antigen presentation to the immune system.
        
        Returns:
        - List of Feature objects representing sampled antigens
        - Empty list if cancer count is below detection threshold has not started sampling before
        """
        antigens = []
        n = self.cancer_n
        
        if n < 300 and not self.antigen_sampling:
            return []  # 300 is the threshold for the immune system to start sampling, after that there's no pausing
        self.antigen_sampling=True
        cancer_grids = self.get_cancer_list()
        
        # Verify cancer count matches list length
        if n != len(cancer_grids):
            print(n, len(cancer_grids), "Cancer count don't match.")
        
        # Sample antigens from random cancer cells
        for i in range(n):
            index = int(random.random() * n)  # Random index
            antigens.append(cancer_grids[index].cell.Feature)
        
        return antigens

# Lymph Node Class

class Lymph_node:
    """
    Represents a lymph node in the immune system, responsible for B cell activation and antibody production.
    
    The lymph node processes antigens from cancer cells, activates matching B cells,
    and produces antibodies that target cancer cells.
    """
    def __init__(self, match_rate, inactive_time, effectiveness):
        """
        Initialize a lymph node with specific parameters.
        
        Parameters:
        - step: Step size for discretizing feature space
        - match_rate: Rate of antigen matching per time step
        - inactive_time: Duration antibodies remain inactive before activation
        - effectiveness: Base effectiveness of produced antibodies
        """
        self.b_cell_reserve = {}  # Dictionary tracking B cell production cooldowns
        self.inactive_time = inactive_time  # Duration antibodies remain inactive
        self.time = 0  # Local time counter
        self.bcell = []  # List of active B cells
        self.escape_cases = []  # Track cases of immune escape
        self.match_rate = match_rate  # Rate of antigen matching
        self.antigens = []  # List of antigens to process
        self.new_anti = 0  # Counter for new antibodies
        self.effectiveness = effectiveness  # Base antibody effectiveness
        
        self.step=0.5 #To allocate B cell features by discretizing the space
    
    def add_antigen(self, antigens):
        """
        Add antigens to the lymph node for processing.
        
        Parameters:
        - antigens: List of Feature objects representing antigens
        """
        self.antigens += antigens
    
    def create_target_immune_cell(self, antigen):
        """
        Create target features for immune cells based on an antigen. Narrow activation
        
        This generates a more focused set of potential B cell features that could
        recognize the antigen, with variations only in one dimension at a time.
        
        Parameters:
        - antigen: Feature object representing the antigen
        
        Returns:
        - List of feature tuples for potential B cells
        """
        # Discretize antigen features to grid points
        optimal_feature = [x - x % self.step for x in antigen.features]
        features = [optimal_feature]
        
        # Create variations by adjusting one dimension at a time
        for i in range(len(optimal_feature)):
            up_feature = [x for x in optimal_feature]
            low_feature = [x for x in optimal_feature]
            up_feature[i] += self.step
            low_feature[i] -= self.step
            features.append(tuple(up_feature))
            features.append(tuple(low_feature))
        
        return features

    def create_broad_target_immune_cell(self, antigen):
        """
        Create a broad range of target features for immune cells based on an antigen.
        
        This generates multiple potential B cell features that could recognize the antigen,
        with variations to cover a broader range.
        
        Parameters:
        - antigen: Feature object representing the antigen
        
        Returns:
        - List of feature tuples for potential B cells
        """
        # Possible neighbor distance
        numbers = [-0.5, 1, 0, -1, 0.5]
        
        # Find a close lower match of B cell with the given antigen
        optimal_feature = [x - x % self.step for x in antigen.features]
        
        # Generate all offset from possible neighbor distance
        combinations = list(itertools.product(numbers, repeat=len(optimal_feature)))
        target = []
        
        # Create all neighbors
        for comb in combinations:
            target.append(tuple([optimal_feature[i] + comb[i] for i in range(len(antigen.features))]))
        
        return target
    def match(self, broad=False):
        """
        Match antigens to B cells and activate them.
        
        This simulates the process of antigen presentation and B cell activation
        in the lymph node.
        
        Parameters:
        - n: Number of B cells to produce per match (default 1)
        - broad: Whether to use broad targeting (default False)
        """
        # Skip if no antigens to process
        if self.antigens == []:
            return
        
        # Determine how many antigens to process in this step
        index = int(self.match_rate * len(self.antigens))
        to_be_matched = self.antigens[:index]  # Take the first portion of antigens
        
        for antigen in to_be_matched:
            # Generate target features based on matching strategy
            if not broad:
                target_immune_cell_features = self.create_target_immune_cell(antigen)
            else:
                target_immune_cell_features = self.create_broad_target_immune_cell(antigen)
            
            
            # Process each target feature
            for feature in target_immune_cell_features:
                feature = tuple(feature)
                try:
                    # Check if B cell is on cooldown
                    if self.b_cell_reserve[feature] == 0:
                        # Reset cooldown and produce B cells
                        self.b_cell_reserve[feature] = 5  # Cooldown period
                        self.bcell.append(Bcell(Feature(feature)))
                except:
                    # Initialize cooldown for new B cell type
                    self.b_cell_reserve[feature] = 0
        
        # Remove processed antigens
        self.antigens = self.antigens[index:]
    
    def diffuse_antibodies(self, Grids, inactive_time, activation):
        """
        Produce antibodies from activated B cells and diffuse them throughout the simulation space.
        
        Parameters:
        - Grids: The Grids object managing the entire simulation space
        - inactive_time: Duration antibodies remain inactive
        - activation: Activation mode affecting antibody multiplier
        """
        # Determine antibody multiplier based on activation mode. For example, instead of activating 3 times more antibodies in TN
        # We simply times a multiplier in the effectiveness function
        n = 1
        if activation == "TN":
            n = 3  # Triple normal effectiveness
        if activation == "DN":
            n = 2  # Double normal effectiveness
        
        grids = Grids.get_list()
        
        # Create antibodies from B cells
        antibodies = [Antibody(Feature(cell.Feature.features), inactive_time, "inactive", self.effectiveness, n) for cell in self.bcell]
        
        # Add antibodies to virtual grid and all real grids
        Grids.virtual_grid.add_antibody(antibodies)
        for grid in grids:
            grid.add_antibody(antibodies)
        
        # Track number of new antibodies
        self.new_anti = len(antibodies) * n
        
        # Clear B cells after antibody production. Basically, once a B cell is matched, its only function is to send one wave of antibody. Then it goes dormant again.
        self.bcell = []  # B cells have short lifespan in this model
    
    def update(self, Grids, activation):
        """
        Update the lymph node for the current time step.
        
        This handles antibody production, antigen processing, and B cell activation/matching.
        
        Parameters:
        - Grids: The Grids object managing the entire simulation space
        - activation: Activation mode affecting antibody production
        """
        # Produce and diffuse antibodies
        self.diffuse_antibodies(Grids,self.inactive_time, activation)
        
        # Sample antigens from cancer cells
        samples = Grids.sample_antigens()
        self.add_antigen(samples)
        
        # Shuffle antigens for random matching
        random.shuffle(self.antigens)
        
        # Match antigens to B cells based on activation mode
        if activation == "B":
            self.match( broad=True)  # Broad matching
        else:
            self.match()  # Normal matching
        
        # Update B cell cooldowns
        for b in self.b_cell_reserve:
            if self.b_cell_reserve[b] > 0:
                self.b_cell_reserve[b] -= 1
        
        # Update local time
        self.time += 1
