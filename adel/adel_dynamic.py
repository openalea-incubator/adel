""" Prototype adel model that uses mtg edition functions"""
from openalea.mtg import MTG
from alinea.adel.AdelR import RunAdel
from openalea.mtg.algo import union
from alinea.adel.astk_interface import AdelWheat, replicate, transform_geom, mtg_interpreter, adel_metamer
from alinea.adel.mtg_editions import find_metamers, add_plant, add_vegetative_metamer, new_mtg_factory
from alinea.adel.AdelR import plantSample


class AdelWheatDyn(AdelWheat):

    def build_stand(self):

        g = MTG()
        plants = self.axeT()
        sample = [int(p) for p in plantSample(self.pars)]
        for i, plant in plants.groupby('plant'):
            plant_properties = {'position': self.positions[i - 1],
                                'azimuth': self.plant_azimuths[i - 1],
                                'refplant_id': sample[i - 1]}
            ms = plant.loc[plant['axe'] == 'MS', :].to_dict('list')
            ms_properties = {'HS_final': float(ms['HS_final'][0]),
                             'nff': int(ms['nf'][0]),
                             'hasEar': bool(int(ms['hasEar'][0])),
                             'azimuth': float(ms['azTb'][0])}
            add_plant(g, i, plant_properties=plant_properties,
                      axis_properties=ms_properties)
        return g

    def add_metamer(self, g, plant=1, axe='MS'):
        vid_plant, vid_axe, metamers = find_metamers(g, plant, axe)
        nff = g.property('nff')[vid_axe]
        num_metamer = len(metamers)
        df = self.phytoT(axe)
        m = df.loc[df['n'] == num_metamer, :].to_dict('list')
        metamer_properties = {'L_shape': m['Ll'][0]}
        ntop = nff - num_metamer + 1
        shape_key = None
        lctype = int(m['Lindex'][0])
        lcindex = int(m['Lseed'][0])
        if lctype != -999 and lcindex != -999:
            shape_key = self.leaves.get_leaf_key(lctype, lcindex, age=None)

        internode_properties = {'ntop': ntop, 'length': 0, 'visible_length': 0,
                                'senesced_length': 0, 'diameter': m['Ed'][0],
                                'azimuth': m['Azim'][0], 'inclination': 0}
        sheath_properties = {'ntop': ntop, 'length': 0, 'visible_length': 0,
                             'senesced_length': 0, 'diameter': m['Gd'][0],
                             'azimuth': 0, 'inclination': 0}
        blade_properties = {'ntop': ntop, 'length': 0, 'visible_length': 0,
                            'rolled_length': 0, 'senesced_length': 0,
                            'n_sect': self.nsect,
                            'shape_mature_length': m['Ll'][0],
                            'shape_max_width': m['Lw'][0],
                            'shape_key': shape_key, 'inclination': 0}

        return add_vegetative_metamer(g, plant, axe, metamer_properties,
                                      internode_properties, sheath_properties,
                                      blade_properties)

    def setup_canopy(self, age=10):

        self.canopy_age = age

        # produce plants positionned at origin
        if self.nrem > 0:
            canopy = RunAdel(age, self.pars_rem, adelpars=self.run_adel_pars)
            grem = new_mtg_factory(canopy, adel_metamer, leaf_sectors=self.nsect,
                               leaves=self.leaves, stand=None, split=self.split,
                               aborting_tiller_reduction=self.aborting_tiller_reduction)
            grem = mtg_interpreter(grem, self.leaves, face_up=self.face_up,
                                   classic=self.classic)

        if self.nquot > 0:
            canopy = RunAdel(age, self.pars_quot, adelpars=self.run_adel_pars)
            gquot = new_mtg_factory(canopy, adel_metamer, leaf_sectors=self.nsect,
                                leaves=self.leaves, stand=None,
                                split=self.split,
                                aborting_tiller_reduction=self.aborting_tiller_reduction)
            gquot = mtg_interpreter(gquot, self.leaves, face_up=self.face_up,
                                    classic=self.classic)
            gquot = replicate(gquot, self.nquot * self.duplicate)

        if self.nrem > 0:
            g = grem
            if self.nquot > 0:
                g = union(g, gquot)
        else:
            g = gquot

        # update positions and domain if smart stand is used
        if self.stand.density_curve is not None:
            new_nplants, self.domain, self.positions, self.domain_area = self.stand.smart_stand(
                self.nplants, at=age)
            self.meta.update(
                {'domain': self.domain, 'domain_area': self.domain_area})
            assert new_nplants == self.nplants

        # dispose plants and renumber them
        pos = g.property('position ')
        az = g.property('azimuth')
        lab = g.property('label')
        geom = g.property('geometry')
        for i, vid in enumerate(g.vertices(1)):
            lab[vid] = 'plant' + str(i + 1)
            pos[vid] = self.positions[i]
            az[vid] = self.plant_azimuths[i]
            for gid in g.components_at_scale(vid, g.max_scale()):
                if gid in geom:
                    geom[gid] = transform_geom(geom[gid], self.positions[i],
                                               self.plant_azimuths[i])
        # add meta
        root = g.node(0)
        self.meta.update({'canopy_age': age})
        root.meta = self.meta
        return g