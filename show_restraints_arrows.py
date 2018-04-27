from pymol import cmd, cgo, CmdException
from pymol import stored
import re
import colorsys

#Usage:
#load prod_0000.pdb; remove solvent; run /home/dimura/opt/pymol_scripts/show_restraints_arrows.py; load_restraints prod_0001.f
def load_restraints(rst_file_path):
  
  DUnames=dict()
  lines = open(rst_file_path, 'r').readlines()
  print(len(lines))
  for il in range (len(lines)):
    line=lines[il]
    #print(il)

    if(line[0]=="#"):
      words = re.split(', | = |;|,| |=',line)
      iat1=words[1].replace("(","").replace(")","")
      iat2=words[4].replace("(","").replace(")","").replace("\n","")
      name1=words[0][1:]
      name2=words[3]
      DUnames[iat1]=name1
      DUnames[iat2]=name2
      #print('# "{}" "{}"'.format(iat1,iat2))
      continue
    
    re1=re.search('iat ?= ?(\d*), ?(\d*),', line)
    if not re1:
      #print('continue not re {}'.format(il))
      continue
    iat1=re1.group(1)
    iat2=re1.group(2)
    r2=re.search('r2 ?= ?(\d*\.?\d*),', line).group(1)
    r3=re.search('r3 ?= ?(\d*\.?\d*),', line).group(1)
    d_target=0.5*(float(r2)+float(r3))
    r4=float(re.search('r4 ?= ?(\d*\.?\d*),', line).group(1))
    r1=float(re.search('r1 ?= ?(\d*\.?\d*),', line).group(1))
    name1=DUnames[iat1]
    name2=DUnames[iat2]
    name=name1+"_"+name2
    if "@" in name:
      name="link_"+name
    else:
      name="rst_"+name
    
    name=name.replace("@","-")
    name=name.replace("'","p")
    name=name.replace("*","s")
    print('{} {} {}'.format(name, iat1, iat2))
    delta_d = cmd.distance(name, "id " + iat1, "id " + iat2) - d_target
    sigma=d_target-r1
    if(delta_d>0):
      sigma=r4-d_target
    deltatos=delta_d/sigma/1.5
    if(abs(deltatos)>1):
      deltatos=deltatos/abs(deltatos)
    hue=0.33+deltatos*0.33
    rgbc=colorsys.hls_to_rgb(hue,0.5,1.0)
    hexc= "0x%0.2X%0.2X%0.2X" % (int(rgbc[0]*255), int(rgbc[1]*255), int(rgbc[2]*255))
    
    cmd.color(hexc,name)
    if(name[:4]=="rst_"):
      arrow_start("id " + iat1,"id " + iat2,5.0*deltatos)
      arrow_start("id " + iat2,"id " + iat1,5.0*deltatos)

      
  cmd.hide("labels", "rst_*")
  cmd.hide("labels", "link_*")
  cmd.set("dash_gap",0.3,"rst_*")
  cmd.set("dash_length",0.3,"rst_*")
  cmd.set("dash_gap",0.1,"link_*")
  cmd.set("dash_length",0.2,"link_*")
  cmd.set("dash_radius",0.04,"link_*")
        
cmd.extend("load_restraints", load_restraints)



def cgo_arrow(atom1='pk1', atom2='pk2', radius=0.5, gap=0.0, hlength=-1, hradius=-1,
              color='blue red', name=''):
    '''
DESCRIPTION
    Create a CGO arrow between two picked atoms.
ARGUMENTS
    atom1 = string: single atom selection or list of 3 floats {default: pk1}
    atom2 = string: single atom selection or list of 3 floats {default: pk2}
    radius = float: arrow radius {default: 0.5}
    gap = float: gap between arrow tips and the two atoms {default: 0.0}
    hlength = float: length of head
    hradius = float: radius of head
    color = string: one or two color names {default: blue red}
    name = string: name of CGO object
    '''
    from chempy import cpv

    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)

    try:
        color1, color2 = color.split()
    except:
        color1 = color2 = color
    color1 = list(cmd.get_color_tuple(color1))
    color2 = list(cmd.get_color_tuple(color2))

    def get_coord(v):
        if not isinstance(v, str):
            return v
        if v.startswith('['):
            return cmd.safe_list_eval(v)
        return cmd.get_atom_coords(v)

    xyz1 = get_coord(atom1)
    xyz2 = get_coord(atom2)
    normal = cpv.normalize(cpv.sub(xyz1, xyz2))

    if hlength < 0:
        hlength = radius * 3.0
    if hradius < 0:
        hradius = hlength * 0.6

    if gap:
        diff = cpv.scale(normal, gap)
        xyz1 = cpv.sub(xyz1, diff)
        xyz2 = cpv.add(xyz2, diff)

    xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

    obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
          [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
          [1.0, 0.0]

    if not name:
        name = cmd.get_unused_name('arrow')

    cmd.load_cgo(obj, name)

cmd.extend('cgo_arrow', cgo_arrow)

def arrow_start(atom1='pk1', atom2='pk2', length=10.0):
    from chempy import cpv
    length=float(length)
    xyz1 = cmd.get_atom_coords(atom1)
    xyz2 = cmd.get_atom_coords(atom2)
    #print xyz1
    #print xyz2
    norm=cpv.normalize(cpv.sub(xyz2, xyz1))
    start=cpv.add(xyz1,cpv.scale(norm,1.4))
    end=cpv.add(cpv.scale(norm,abs(length)),start)
    if(length<0.0):
      cgo_arrow(end,start,radius=0.3,color='red')
    else:
      cgo_arrow(start,end,radius=0.3,color='blue')
    
cmd.extend('arrow_start', arrow_start)