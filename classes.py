#!/usr/bin/env python3
import numpy as np
import pandas as pd
import netCDF4 as nc


class Fort14():
    def __init__(self, fname):
        """Initialize from a given path to fort.14
        """
        df = pd.read_csv(fname, sep="\s+", names=list('abcde'),on_bad_lines='skip')
        # Convert to appropriate types
        df[['b', 'c', 'd', 'e']] = df[['b', 'c', 'd', 'e']].astype('object')
        df['b'][1] = int(df['b'][1])
        df['a'][1] = int(df['a'][1])

        self.df = df
        self.num_nodes = self.df['b'][1]
        self.num_elems = self.df['a'][1]

        self.df['a'][2:2+self.num_nodes] = self.df['a'][2:2+self.num_nodes].astype(int)
        self.node_indices = self.df['a'][2:2+self.num_nodes].to_numpy()

        self.x = self.df['b'][2:2+self.num_nodes].to_numpy()
        self.y = self.df['c'][2:2+self.num_nodes].to_numpy()

        # Element connectivity
        offset = 2 + self.num_nodes
        self.df.iloc[offset:offset+self.num_elems, 1:] = self.df.iloc[offset:offset+self.num_elems, 1:].astype(int)
        self.connectivity = self.df.iloc[offset:offset+self.num_elems,2:].to_numpy()

        # Get true node numbering
        self.true_node_indices, self.true_connectivity = self.renumber()

        # TODO: renumber boundary nodes

    def bounding_box(self):
        return [np.min(self.x), np.max(self.x), np.min(self.y), np.max(self.y)]

    def is_in_mesh(self, x, y):
        """Return True if the point x,y is contained in the mesh.
        """
        min_x = np.min(self.x)
        max_x = np.max(self.x)

        min_y = np.min(self.y)
        max_y = np.max(self.y)

        return (x >= min_x and x <= max_x) and (y >= min_y and y <= max_y)

    def new_index(self, node):
        """Return the true index of a given node number."""
        return np.where(self.node_indices == node)[0][0] + 1

    def renumber(self):
        """Return nodes renumbered starting from 1."""
        new_index_v = np.vectorize(self.new_index)
        true_node_indices = new_index_v(self.node_indices)
        true_connectivity = np.array([new_index_v(tri) for tri in self.connectivity])

        assert np.min(true_node_indices) == np.min(true_connectivity[:])
        assert np.max(true_node_indices) == np.max(true_connectivity[:])

        return true_node_indices, true_connectivity

    def true_numbering(self):
        """Return a mesh DataFrame with true node numbering."""
        copy = self.df.copy()
        copy['a'][2:2+self.num_nodes] = self.true_node_indices

        offset = 2 + self.num_nodes

        #copy.iloc[offset:offset+self.num_elems, 1] = copy.iloc[offset:offset+self.num_elems, 1].astype(int)
        copy.iloc[offset:offset+self.num_elems, 2:5] = self.true_connectivity

        return copy

    def write(self, fname):
        """Write the true numbering to a new file."""
        copy = self.true_numbering()
        copy.to_csv(fname, index=False, header=False, sep="\t")


#f14 = Fort14("120m_nolevee/PE0571/fort.14")
