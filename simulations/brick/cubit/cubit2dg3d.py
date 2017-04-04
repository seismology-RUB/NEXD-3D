import cubit
import numpy
from math import sqrt

class mesh():
    def __init__(self):
        print('pass')
        # set main variables
        self.mesh_name='mesh'
        self.coord_name='coord'
        self.tet_name='tet'
        self.face_name='face'
        self.neighbor_name='neighbor'
        self.material_name='mat'
        self.matprop_name='matprop'
        self.absorb_name='absorb'
        self.free_name='free'
        self.tet='TETRA4'
        self.tri='TRISHELL3'

        self.get_mesh()

        self.get_block()
#        self.get_face()

        self.coord_write(self.coord_name)
        self.elem_write(self.mesh_name)
        self.mat_write(self.material_name)
        self.matprop_write(self.matprop_name)
        self.absorb_write(self.absorb_name)
        self.free_write(self.free_name)


    def get_block(self):
        ''' extract block information '''
        block_flag=[]
        block_mat=[]
        block_free=[]
        block_absorb=[]
        block_pml=[]
        material=[]
        # number of blocks
        blocks=cubit.get_block_id_list()

        # go through blocks
        for block in blocks:
            name=cubit.get_exodus_entity_name('block',block)
            type=cubit.get_block_element_type(block)
            print block,name,blocks,type

            if type == self.tet:
                # check if elastic
                if name.find('elastic') >=0:
                    print('found elastic block')
                    imat=1
                else:
                    print('error, no correct material in block',block)
                
                # get number of attributes
                nattrib=cubit.get_block_attribute_count(block)
                if nattrib != 0:
                    flag=int(cubit.get_block_attribute_value(block,0))
                    vp=cubit.get_block_attribute_value(block,1)
                    vs=cubit.get_block_attribute_value(block,2)
                    rho=cubit.get_block_attribute_value(block,3)
                    qp=cubit.get_block_attribute_value(block,4)
                    qs=cubit.get_block_attribute_value(block,5)
                    pml=int(cubit.get_block_attribute_value(block,6))
                else:
                    print('error no attributes in block',block)
                block_flag.append(flag)
                block_mat.append(block)
                block_pml.append(pml)
                material.append([imat,vp,vs,rho,qp,qs])
            elif type == self.tri :
                if name.find('free') >=0:
                    print('found free surface')
                    block_free.append(block)
                elif name.find('absorb') >=0:
                    print('found absorb surface')
                    block_absorb.append(block)
                else:
                    print('error in boundaries',block)

#        ns=cubit.get_nodeset_id_list()
#        for n in ns:
#            name=cubit.get_exodus_entity_name('nodeset',n)
#            if name.find('free') >=0:
#                print('found free surface nodes')
#                ns_free.append(n)
#            elif name.find('absorb') >=0:
#                print('found absorb surface nodes')
#                ns_absorb.append(n)
#            else:
#                print('error in boundaries',n)
            
        print('BLOCKFLAG',block_flag)
        self.block_flag=block_flag
        self.block_mat=block_mat
        self.mat=material
        self.block_free=block_free
        self.block_absorb=block_absorb
        self.block_pml=block_pml


    def get_mesh(self):
        ''' get tri mesh from cubit in format for dg3d '''
        print('Start extracting mesh for dg3d')

    def coord_write(self,coord_name):
        ''' write nodes file '''
        coord=open(coord_name,'w')
        print('writing '+coord_name)
        node_list=cubit.parse_cubit_list('node','all')
        ncoord=len(node_list)
        print('number of nodes:',str(ncoord))
        # write ncoord
        coord.write('%10i\n' % ncoord)
        
        #write coords
        for node in node_list:
            x,y,z=cubit.get_nodal_coordinates(node)
            txt=('%10i %20f %20f %20f\n') % (node,x,y,z)
            coord.write(txt)
        coord.close()
        print('Ok')

    def mat_write(self,mat_name):
        mat=open(mat_name,'w')
        print('Writing '+mat_name+'.....')
        for block,flag,pml in zip(self.block_mat,self.block_flag,self.block_pml):
                tets=cubit.get_block_tets(block)
                for tet in tets:
                    mat.write(('%10i %10i %10i\n') % (tet,flag,pml))
        mat.close()
        print('Ok')

    def elem_write(self,mesh_name):
        meshfile=open(mesh_name,'w')
        print('Writing '+mesh_name+'.....')
        nelem=cubit.get_tet_count()
        print('number of elements:',str(nelem))
        meshfile.write(str(nelem)+'\n')
        num_write=0
        temp_tet=[]
        tetlength=1
        for block,flag in zip(self.block_mat,self.block_flag):
            tetlength += len(cubit.get_block_tets(block))
        tet_vp=range(tetlength)
        tet_vs=range(tetlength)
        tet_rho=range(tetlength)
        tet_qp=range(tetlength)
        tet_qs=range(tetlength)
        tet_block=range(tetlength)
        for block,flag in zip(self.block_mat,self.block_flag):
            tets=cubit.get_block_tets(block)
            for tet in tets:
                nodes=cubit.get_connectivity('Tet',tet)
                tet_vp[tet]=cubit.get_block_attribute_value(block,1)
                tet_vs[tet]=cubit.get_block_attribute_value(block,2)
                tet_rho[tet]=cubit.get_block_attribute_value(block,3)
                tet_qp[tet]=cubit.get_block_attribute_value(block,4)
                tet_qs[tet]=cubit.get_block_attribute_value(block,5)
                tet_block[tet]=block
                temp_tet.append(tet)
#                txt=('%10i ')% tet
#                txt=txt+('%10i %10i %10i %10i\n')% nodes[:]
#                meshfile.write(txt)

        temp_tet.sort()
        for tet in temp_tet:
            nodes=cubit.get_connectivity('Tet',tet)
            txt=('%10i ')% tet
            txt=txt+('%10i %10i %10i %10i') % nodes[:]
            txt=txt+(' %9.1f %9.1f %9.1f %5i %5i\n') % (tet_vp[tet],tet_vs[tet],tet_rho[tet],tet_qp[tet],tet_qs[tet])
            meshfile.write(txt)

        meshfile.close()
        print('Ok')

    def matprop_write(self,matprop_name):
        matpropfile=open(matprop_name,'w')
        print('writing '+matprop_name)
        nmat=len(self.block_mat)
        print('number of materials: ',str(nmat))
        matpropfile.write(str(nmat)+'\n')
        for i,block in enumerate(self.block_mat):
#            print(block)
#            print(block,self.mat[i][:])
            txt=str(block)
            txt=txt+' '+str(self.mat[i][0])+' '+str(self.mat[i][1])+' '+str(self.mat[i][2])+' '+str(self.mat[i][3])+' '+str(self.mat[i][4])+' '+str(self.mat[i][5])
            matpropfile.write(txt+'\n')
        
#    def get_face(self):
#        for n in self.block_free:
#            nodes_free=cubit.get_nodeset_nodes(n)
#            print(nodes_free)

    def absorb_write(self,absorb_name):
        absorbfile=open(absorb_name,'w')
        print('writing '+absorb_name)
        
#        for n in self.ns_absorb:
#            nodes_absorb=cubit.get_nodeset_nodes(n)
#        nabs=len(nodes_absorb)
#        absorbfile.write(str(nabs)+'\n')
#        print('Number of absorbing nodes ',nabs)

#        for node in nodes_absorb:
#            absorbfile.write(str(node)+'\n')
###
        list_tet=cubit.parse_cubit_list('tet','all')
        
        if (self.block_absorb != []):
	    for n in self.block_absorb:
                elems_absorb=cubit.get_block_tris(n)
                dic_absorb=dict(zip(elems_absorb,elems_absorb))
                nabs=len(elems_absorb)
                absorbfile.write(str(nabs)+'\n')
                print('Number of absorbing elements ',nabs)

            for tet in list_tet:
                tris=cubit.get_sub_elements('tet',tet,2)
                for tri in tris:
                    if dic_absorb.has_key(tri):
                        nodes=cubit.get_connectivity('Tri',tri)
                        txt=('%10i ')% tet
                        txt=txt+('%10i %10i %10i')% nodes[:]
                        absorbfile.write(str(txt)+'\n')
	else:
	    absorbfile.write(str(0)+'\n')
            print('Number of absorbing elements ',0)

#        for tri in elems_absorb:
#            nodes=cubit.get_connectivity('Tri',tri)
#            txt=('%10i ')% tri
#            txt=txt+('%10i %10i %10i')% nodes[:]
#            absorbfile.write(str(txt)+'\n')

###

    def free_write(self,free_name):
        freefile=open(free_name,'w')
        print('writing '+free_name)

        list_tet=cubit.parse_cubit_list('tet','all')
        
        if (self.block_free != []):
            for n in self.block_free:
                elems_free=cubit.get_block_tris(n)
                dic_free=dict(zip(elems_free,elems_free))
                nfree=len(elems_free)
                freefile.write(str(nfree)+'\n')
                print('Number of free elements ',nfree)
        
            for tet in list_tet:
                tris=cubit.get_sub_elements('tet',tet,2)
                for tri in tris:
                    if dic_free.has_key(tri):
                        nodes=cubit.get_connectivity('Tri',tri)
                        txt=('%10i ')% tet
                        txt=txt+('%10i %10i %10i')% nodes[:]
                        freefile.write(str(txt)+'\n')
        else:
            freefile.write(str(0)+'\n')
            print('Number of free elements ',0)


if __name__ == '__main__':
    mesh()
