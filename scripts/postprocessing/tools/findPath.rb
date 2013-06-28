#!/usr/bin/env ruby
require 'pp'
require 'facets/array/nonuniq'
require 'numru/netcdf'
include NumRu

module IconPathsAlongCells
  # R2B04: list of cellsof which their common edge define the path
  CellPairsLists= {
    :gibraltar      => [[4543,4485]],

    :denmarkStrait  => [[4690,4695]],

    :drakePassage   => [
      [18481,18495],
      [18483,18493],
      [18642,18649],
      [18651,18643],
      [19186,19191]],

    :indonesianThroughflow_00 => [
      [6864,6866],
      [6868,6869],
      [6872,6875],
      [12064,12067],
      [12072,12073],
      [12152,12154]],

    :indonesianThroughflow_01 => [
      [6855,6857],
      [6858,6869],
      [6878,6879],
      [11265,11269],
      [11461,11462],
      [11470,11471],
      [11509,11510],
      [11518,11519],
      [12106,12107],
      [12109,12110],
      [12133,12134],
      [12142,12143]],

    :scotland_iceland => [
      [4238,4256],
      [4259,4268],
      [4270,4271],
      [4682,4702]],
  }

  def IconPathsAlongCells.getEdgesAndVerts(ifile)
    ifileHandle  = NetCDF.open(ifile,"r")
    cellEdges    = ifileHandle.var('edge_of_cell').get
    cellVertices = ifileHandle.var('vertex_of_cell').get
    edgeVertices = ifileHandle.var('edge_vertices').get

    paths = {}

    puts '#==========================================================================='
    CellPairsLists.each {|location,cellPairs|
      next unless :scotland_iceland == location
      puts '============================================================================='
      puts location
      puts ["cells: ", cellPairs.join(' ')].join

      verts, edges = [], []
      cellPairs.each {|cellpair|
          verts = []
        cellpair.each {|cell|
          edges << cellEdges[cell,0..-1].to_a
          verts << cellVertices[cell,0..-1].to_a
        }
      }
      verts = verts.flatten.nonuniq
      edges = edges.flatten.nonuniq

      paths[location] = {verts: verts, edges: edges}

      puts ["common verts: ", verts.join(',')].join
      puts ["common edges: ", edges.join(',')].join
    }
    puts '============================================================================='

    return paths
  end
end

if $0 == __FILE__ then
  # = MAIN =======================================================================
  ifile = ARGV[0]

  if ifile.nil? or not File.exist?(ifile) 
    warn "Cound not read input file '#{ifile}'!"
    warn "Usage:\n\t./findPath.rb <icon-grid-file>"
    exit(1)
  end

  a = IconPathsAlongCells.getEdgesAndVerts(ifile)
end
