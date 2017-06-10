"""
Usage:
genome_visualization.py (-view FILE) [(-ref NAME --start VALUE --end VALUE )]  [-reff FASTAFILE] [--mp]
genome_visualization.py (-h | --help)

Visualizes or plots reads in SAM/BAM files to aid in the analysis of data through different
plotting manners with features enhancing analyzing effects


Arguments:
  FILE      input file to be visualized [BAM format]
  VALUE     start and end regions' values
  <path>      directory at which to execute the plots
  FASTAFILE   input fasta refeerence genome


Options:
  -view FILE                   Imports and Views BAM files
  -ref  NAME                   Name of the reference chromosome
  -s,--start VALUE             Determines specific starting region in file
  -e, --end VALUE              Determines the ending region in file (optional)[default: end]
  --mp                         Displays matepairs among plots
  -h, --help                   Shows help document and quit
  -reff FASTAFILE              Plots reference genome

"""

from docopt import docopt
import numpy as np
import pysam
import sys
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
from Bio import SeqIO
from plotly import tools

def func(zz,ref_ann,value):
    if value=='C':
        zz.append(0.1)
        ref_ann.append('C')
    elif value=='G':
        zz.append(0.2)
        ref_ann.append('G')
    elif value=='A':
        zz.append(0.3)
        ref_ann.append('A')
    elif value=='T':
        zz.append(0.4)
        ref_ann.append('T')
    return (zz, ref_ann)


def pairedend(traces,  x1,x2,y1,y2):

    end=len(x1) +2
    z1= zip(x1,x2)
    z2=zip(y1,y2)
    counter=0

    x_data = [
                [0,0, 0, 0],
                [0,0, 0, 0],
                [0,0, 0, 0],
    ]

    y_data = [
                [0,0, 0, 0],
                [0,0, 0, 0],
                [0,0, 0, 0],
    ]


    for x1,x2 in z1:
        for y1,y2 in z2:
           y_data.append([0, counter, counter, 0])
           counter=counter-0.2
           if counter < -3:
               counter =0
        x_data.append([x1,x1,x2,x2])


    for i in range(3, end):
        traces.append(go.Scatter(
            x=x_data[i],
            y=y_data[i],
            mode='lines',
            name="'paired match",
            line=dict(color=('rgb(205, 12, 24)'),
                      dash='dot'),
            showlegend=False,
            connectgaps=True,
        ))

    return traces

def ref_plot(traces,reffile, start, end):
    old=0
    new=0
    zz=[]
    int(start)
    int(end)
    reflength=end- start
    refcounter=0
    ref_ann=[]

    try:
        input_file = open(reffile, 'r')
        for cur_record in SeqIO.parse(input_file, "fasta") :
            gene_name = cur_record.name
            for item in cur_record[start:]:
                if refcounter ==reflength:
                    break
                if refcounter==0:
                    refcounter+=1
                    continue
                if old==0:
                    old=item
                    zz, ref_ann= func(zz,ref_ann,item)
                    refcounter+=1
                    continue
                new= item
                zz, ref_ann= func(zz,ref_ann, new)
                old=new
                refcounter+=1
        input_file.close()

        data=[]

        x= np.arange(start+1, end+1,1)
        data.append(zz)
        #colorscale=[[0.1, '#ff0000'], [.2, '#00ff00'],[.3, '#0000ff'],[0.4, '#000000']]
        colorscale="Jet"
        traces.append(go.Heatmap(z=data,x=x, showlegend=False,text=ref_ann, colorscale=colorscale,showscale=False,  zmin=0.1, zmax=0.4))
        return (gene_name, traces)
    except:
        return ("No Reference Genome Entered", traces)


def plot(file, reference, start, stop,refgenome):
  end=[]
  position=[]
  levels=[-5]
  found=False
  x1=[]
  y1=[]
  x2=[]
  y2=[]
  i=0
  text_rev=[]
  text_forw=[]
  number_of_reads=[]
  x_rev=[[0,0],
        [0,0]]
  y_rev=[[0,0],
        [0,0]]
  x_forw=[[0,0],
        [0,0]]
  y_forw=[[0,0],
        [0,0]]
  traces=[]

  start= int(start)
  stop= int(stop)
  max_coeff=0
  final_coeff=0


  for read in file.fetch(reference,start , stop):
      end.append(read.reference_end)
      position.append(read.pos)
      if i==0: position.pop()

#create annotations text
      if read.is_paired ==True: pair_status='Read is paired.'
      else: pair_status='Read not paired.'

      if read.is_qcfail ==True: QC_status='Failed QC: Yes'
      else: QC_status='Failed QC: No'

      if read.is_secondary==True: alignment_status='Not Primary Alignment'
      else: alignment_status='Primary Alignment'
      tags=[]
      for tag in read.tags:
          tags.append(tag)


##continue with plot function
      if read.is_paired:
          if read.is_proper_pair:
              pos1= read.pos
              #name= read.query_name
              end1=read.reference_end
              point1=(pos1+ end1) /2
              mate= file.mate(read)
              pos2= mate.pos
              end2=mate.reference_end
              point2=(pos2 + end2 )/2
              if i>0:
                  x1.append(point1)
                  y1.append(0)
                  x2.append(point2)
                  y2.append(0)


      for m,l in enumerate(levels):

          if read.pos > l:
              i = read.reference_end
              levels[m]=i
              found=True
              coeff=1
              break

          if not found:
              i= read.reference_end
              levels.append(i)

          x=read.pos
          yaxis= coeff
          if read.flag==16:
              x_rev.append([x,i])
              y_rev.append([yaxis,yaxis])
              text_rev.append('\nRead Name: %s \nRead Length: %s: \nRead Sequence: %s \nRead Flag: %s \nReference ID: %s \nReference Start: %s '
                  '\nNext Reference ID: %s \nNext Reference Start: %s \nMapping Quality: %s \nCIGAR: %s \nTemplate Length: %s \n%s \n%s \n%s \nTags: \n%s' % (
                      read.query_name, read.query_alignment_length, read.query_sequence,  read.flag, read.reference_id, read.reference_start, read.next_reference_id,
                      read.next_reference_start, read.mapping_quality, read.cigar, read.template_length, pair_status, QC_status, alignment_status, tags))

          else:
              x_forw.append([x,i])
              y_forw.append([yaxis,yaxis])
              text_forw.append('\nRead Name: %s \nRead Length: %s: \nRead Sequence: %s \nRead Flag: %s \nReference ID: %s \nReference Start: %s '
                  '\nNext Reference ID: %s \nNext Reference Start: %s \nMapping Quality: %s \nCIGAR: %s \nTemplate Length: %s \n%s \n%s \n%s \nTags: \n%s' % (
                      read.query_name, read.query_alignment_length, read.query_sequence,  read.flag, read.reference_id, read.reference_start, read.next_reference_id,
                      read.next_reference_start, read.mapping_quality, read.cigar, read.template_length, pair_status, QC_status, alignment_status, tags))


          coeff=coeff+ 1
          final_coeff= coeff
          if final_coeff>max_coeff:
              max_coeff =final_coeff

          if read.is_paired:
              if read.is_proper_pair:
                  y1.pop()
                  y1.append(yaxis)


          for index,value in enumerate(x2):
              if x==value:
                  y2[index]=yaxis
  max_numb=0
  final_numb=0
  for pileupcolumn in file.pileup(reference, start, stop):
      numofreads=pileupcolumn.n
      number_of_reads.append(numofreads)

      final_num= numofreads
      if final_numb>max_numb:
          max_numb =final_numb


  length_rev=len(x_rev)
  length_forw=len(x_forw)
  traces.append(go.Scatter(
      x=x_rev,
      y=y_rev,
                mode='lines',
                name="reverse read",
        textposition="topright",
                textfont=dict(family='Arial', size=5, color='#ffffff'),
                                line=dict(
                    color=('#ff7f50'),
                    width=4,)
            ))

  traces.append(go.Scatter(
                  x=x_forw,
                  y=y_forw,
                  mode='lines',
                  name="forward read",
        textposition="topright",
                  textfont=dict(family='Arial', size=5, color='#ffffff'),

                  line=dict(
                  color=('#10160a'),
                  width=4,
                  )
              ))

  for mmn in range(2, length_rev):
    traces.append(go.Scatter(
                x=x_rev[mmn],
                y=y_rev[mmn],
                mode='lines',
                     textposition="topright",
                textfont=dict(family='Arial', size=5, color='#ffffff'),
                text=text_rev[mmn-2],
                line=dict(
                    color=('#ff7f50'),
                    width=4,),
        showlegend=False
            ))
  for mmm in range(2, length_forw):
    traces.append(go.Scatter(
                  x=x_forw[mmm],
                  y=y_forw[mmm],
                  mode='lines',
        textposition="topright",
                  textfont=dict(family='Arial', size=5, color='#ffffff'),
                  text=text_forw[mmm-2],
                  line=dict(
                  color=('#10160a'),
                  width=4,),
        showlegend=False
              ))


  traces.append(go.Bar(
          x = position ,
          y = number_of_reads,
      yaxis='y2',
      xaxis='x2',
      name='Coverage',
      marker=dict(
        color='rgb(49,130,189)'
    )))

  gene_name, traces= ref_plot(traces,refgenome, start, i)


  layout = go.Layout(
                font=dict(family='Courier New, monospace', size=15, color='#7f7f7f'),
                titlefont=dict(
                    family='Courier New, monospace',
                    size=30,
                    color='#7f7f7f'
                    ),
                     title = 'Genome Visualization ',

                     paper_bgcolor='rgb(255,255,255)',

                     plot_bgcolor='rgb(229,229,229)',

              xaxis=dict(

                  title='Genome Sequence',
                  titlefont=dict(
                  family='Courier New, monospace',
                  size=18,
                  color='#7f7f7f'
                   ),
              gridcolor='rgb(255,255,255)',
              showgrid=True,
        #      showline=False,
              showticklabels=True,
              tickcolor='rgb(127,127,127)',
              ticks='outside',
              zeroline=False
              ),

      yaxis=dict(range=[-3.5,max_coeff+2],
        domain=[0, 0.85]
    ),
    legend=dict(
        traceorder='reversed'
    ),
    yaxis2=dict(
        domain=[0.85, 1],
        zeroline=False,
        ticktext = [
                "Read Coverage"
            ],
        tickvals=[max_numb]
  ),
      xaxis2=dict(
        anchor='y2',
          #showticklabels=False
    ),
           height=800,
          width=1000,
          autosize=False,

          margin=dict(
              autoexpand=False,
              l=180,
              r=-5,
              t=100,
              b=150,
              pad=7
          ),

      )

  annotations=[]

  annotations.append(dict(
            x=start+((i-start)/2),
            y=-3,
            xanchor='center',
            yanchor='top',
            text='Genome: %s' % (gene_name),
            font=dict(family='Arial',
            size=12,
            color='rgb(150,150,150)'),
            showarrow=False))

  y_label=["Ref Genome", "Read Alignment","Read Coverage", "Tracks"]
  y_trace=[ 0, max_coeff/2, max_coeff+5, max_coeff+final_coeff+5]
  for y_label,y_trace in zip(y_label, y_trace):
      annotations.append(
      dict(
      xref='paper',
          x=-0.03,
          y=y_trace,
          xanchor='right',
          yanchor='middle',
          text=y_label,
          font=dict(family='Arial',
                    size=18),
          showarrow=False))


  layout['annotations'] = annotations

  if (pairs): data=pairedend (traces,x1, x2, y1, y2)
  else: data=traces

  fig = go.Figure(data=data, layout=layout)
  plotly.offline.plot(fig)


def plotall(file, refgenome):
  end=[]
  position=[]
  levels=[-5]
  found=False
  x1=[]
  y1=[]
  x2=[]
  y2=[]
  i=0
  #pos=[]
  number_of_reads=[]
  text_rev=[]
  text_forw=[]

  x_rev=[[0,0],
        [0,0]]
  y_rev=[[0,0],
        [0,0]]
  x_forw=[[0,0],
        [0,0]]
  y_forw=[[0,0],
        [0,0]]
  traces =[]

  for read in file.fetch(until_eof=True):
      end.append(read.reference_end)
      position.append(read.pos)
      if i==0: position.pop()


#create annotations text
      if read.is_paired ==True: pair_status='Read is paired.'
      else: pair_status='Read not paired.'

      if read.is_qcfail ==True: QC_status='Failed QC: Yes'
      else: QC_status='Failed QC: No'

      if read.is_secondary==True: alignment_status='Not Primary Alignment'
      else: alignment_status='Primary Alignment'
      tags=[]
      for tag in read.tags:
          tags.append(tag)


##continue with plot function

      if read.is_paired:
          if read.is_proper_pair:
              pos1= read.pos
              end1=read.reference_end
              point1=(pos1+ end1) /2

              mate= file.mate(read)
              pos2= mate.pos
              end2=mate.reference_end
              point2=(pos2 + end2 )/2
              if i>0:
                  x1.append(point1)
                  y1.append(0)
                  x2.append(point2)
                  y2.append(0)


      for m,l in enumerate(levels):
          if read.pos > l:
              i = read.reference_end
              levels[m]=i
              found=True
              coeff=1
              break

          if not found:
              i= read.reference_end
              levels.append(i)

          x=read.pos
          yaxis= coeff
          if read.flag==16:
              x_rev.append([x,i])
              y_rev.append([yaxis,yaxis])
              text_rev.append('\nRead Name: %s \nRead Length: %s: \nRead Sequence: %s \nRead Flag: %s \nReference ID: %s \nReference Start: %s '
                  '\nNext Reference ID: %s \nNext Reference Start: %s \nMapping Quality: %s \nCIGAR: %s \nTemplate Length: %s \n%s \n%s \n%s \nTags: \n%s' % (
                      read.query_name, read.query_alignment_length, read.query_sequence,  read.flag, read.reference_id, read.reference_start, read.next_reference_id,
                      read.next_reference_start, read.mapping_quality, read.cigar, read.template_length, pair_status, QC_status, alignment_status, tags))

          else:
              x_forw.append([x,i])
              y_forw.append([yaxis,yaxis])
              text_forw.append('\nRead Name: %s \nRead Length: %s: \nRead Sequence: %s \nRead Flag: %s \nReference ID: %s \nReference Start: %s '
                  '\nNext Reference ID: %s \nNext Reference Start: %s \nMapping Quality: %s \nCIGAR: %s \nTemplate Length: %s \n%s \n%s \n%s \nTags: \n%s' % (
                      read.query_name, read.query_alignment_length, read.query_sequence,  read.flag, read.reference_id, read.reference_start, read.next_reference_id,
                      read.next_reference_start, read.mapping_quality, read.cigar, read.template_length, pair_status, QC_status, alignment_status, tags))


          coeff=coeff+ 1
          final_coeff= coeff
          if final_coeff>max_coeff:
              max_coeff =final_coeff

          if read.is_paired:
              if read.is_proper_pair:
                  y1.pop()
                  y1.append(yaxis)


          for index,value in enumerate(x2):
              if x==value:
                  y2[index]=yaxis

  max_numb=0
  final_numb=0
  for pileupcolumn in file.pileup(until_eof=True):
      numofreads=pileupcolumn.n
      number_of_reads.append(numofreads)

      final_num= numofreads
      if final_numb>max_numb:
          max_numb =final_numb


  length_rev=len(x_rev)
  length_forw=len(x_forw)
  traces.append(go.Scatter(
      x=x_rev,
      y=y_rev,
                mode='lines',
                name="reverse read",
        textposition="topright",
                textfont=dict(family='Arial', size=5, color='#ffffff'),
                                line=dict(
                    color=('#ff7f50'),
                    width=4,)
            ))

  traces.append(go.Scatter(
                  x=x_forw,
                  y=y_forw,
                  mode='lines',
                  name="forward read",
        textposition="topright",
                  textfont=dict(family='Arial', size=5, color='#ffffff'),

                  line=dict(
                  color=('#10160a'),
                  width=4,
                  )
              ))


  for mmn in range(2, length_rev):
    traces.append(go.Scatter(
                x=x_rev[mmn],
                y=y_rev[mmn],
                mode='lines',
                     textposition="topright",
                textfont=dict(family='Arial', size=5, color='#ffffff'),
                text=text_rev[mmn-2],
                line=dict(
                    color=('#ff7f50'),
                    width=4,),
        showlegend=False
            ))
  for mmm in range(2, length_forw):
    traces.append(go.Scatter(
                  x=x_forw[mmm],
                  y=y_forw[mmm],
                  mode='lines',
        textposition="topright",
                  textfont=dict(family='Arial', size=5, color='#ffffff'),
                  text=text_forw[mmm-2],
                  line=dict(
                  color=('#10160a'),
                  width=4,),
        showlegend=False
              ))

       # Create a trace for the positions plot
  traces.append(go.Bar(
          x = position ,
          y = number_of_reads,
      yaxis='y2',
      xaxis='x2',
      name='Coverage',
      marker=dict(
        color='rgb(49,130,189)'
    )))

  gene_name, traces= ref_plot(traces,refgenome, start, i)




  layout = go.Layout(font=dict(family='Courier New, monospace', size=15, color='#7f7f7f'),
                titlefont=dict(
                    family='Courier New, monospace',
                    size=30,
                    color='#7f7f7f'
                    ),
                     title = 'Genome Visualization ',

                     paper_bgcolor='rgb(255,255,255)',

                     plot_bgcolor='rgb(229,229,229)',

                     xaxis=dict(

                  title='Genome Sequence',
                  titlefont=dict(
                  family='Courier New, monospace',
                  size=18,
                  color='#7f7f7f'
                   ),
              gridcolor='rgb(255,255,255)',
              showgrid=True,
        #      showline=False,
              showticklabels=True,
              tickcolor='rgb(127,127,127)',
              ticks='outside',
              zeroline=False
              ),
      yaxis=dict(range=[-3.5,max_coeff+2],
        domain=[0, 0.85]
    ),
    legend=dict(
        traceorder='reversed'
    ),
    yaxis2=dict(
        domain=[0.85, 1],
        zeroline=False,
        ticktext = [
                "Read Coverage"
            ],
        tickvals=[max_numb]
  ),
      xaxis2=dict(
        anchor='y2',
          #showticklabels=False
    ),
           height=800,
          width=1000,
          autosize=False,

          margin=dict(
              autoexpand=False,
              l=180,
              r=-5,
              t=100,
              b=150,
              pad=7
          ),

      )
  annotations=[]

  annotations.append(dict(
            x=start+((i-start)/2),
            y=-3.25,
            xanchor='center',
            yanchor='top',
            text='Genome: %s' % (gene_name),
            font=dict(family='Arial',
            size=12,
            color='rgb(150,150,150)'),
            showarrow=False))

  y_label=["Ref Genome", "Read Alignment","Read Coverage", "Tracks"]
  y_trace=[ 0, max_coeff/2, max_coeff+5, max_coeff+final_coeff+5]
  for y_label,y_trace in zip(y_label, y_trace):
      annotations.append(
      dict(
      xref='paper',
          x=-0.03,
          y=y_trace,
          xanchor='right',
          yanchor='middle',
          text=y_label,
          font=dict(family='Arial',
                    size=18),
          showarrow=False))


  layout['annotations'] = annotations

  if (pairs): data=pairedend (traces,x1, x2, y1, y2)
  else: data=traces

  fig = go.Figure(data=data, layout=layout)
  plotly.offline.plot(fig)


if __name__ == '__main__':
    arguments = docopt(__doc__)

    filename = arguments['FILE']
    file= pysam.AlignmentFile(filename)

    refgenome = arguments['FASTAFILE']
    pairs = arguments['--mp']
    start = arguments['--start']
    stop = arguments['--end']
    reference=arguments['NAME']
    #path=arguments['<path>']


    if (start):
        stop=stop[2]
        plot(file, reference, start, stop,refgenome)

    else:
        plotall(file,refgenome)
    file.close()