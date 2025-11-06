from __future__ import annotations
import argparse
import re
import sys
from typing import Dict, List, Tuple
from collections import defaultdict

#Tuple for coordinate
Coord = Tuple[float, float, float]


def load_points(path: str):
    #Load, points, read 'x y z' into {id: (x,y,z)}
    pts: Dict[int, Coord] = {}
    with open(path, "r", encoding="utf-8") as f:
        #Read line and split
        for ln, raw in enumerate(f, 0):
            s = raw.strip()

            parts = s.split()
            
            #Not enough points
            if len(parts) < 3:
                # Skip line
                continue
            
            #Convert to coordinates x, y, z
            try:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])

            #Skip line that doesn't parse cleanly
            except ValueError:

                continue
            
            #Add to the dictionary
            pts[ln] = (x, y, z)
    return pts


def load_clusters(path: str) -> Dict[int, List[int]]:
    #Parse clusters file into {centroid_id: [point_ids...]}.
    #Format:
    #    7:
    #    81 19 27

    #Each centroid may store multiple values, use default dictionary
    clusters = defaultdict(list)
    current = None

    #Read file
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            
            #String starts with #, skip
            if not s or s.startswith("#"):
                continue

            #Char : found in s
            if ":" in s:
                #Split string
                left, right = s.split(":", 1)
                left = left.strip()
                
                #Set current centroid id
                try:
                    current = int(left)
                
                #Impossible, throw an exception
                except ValueError:
                    current = None
                    continue

            #Next line after the char?
            else:
                #Line is empty
                if current is None:
                    continue
                
                #Line contains clients
                for tok in s.split():
                    
                    #Process one by one
                    try:
                        clusters[current].append(int(tok))
                    except ValueError:
                        pass

    return dict(clusters)

def load_centroid_ids(path: str, col:int):
    #Read centroid ids from the given column
    ids: List[int] = []
    #Open file
    with open(path, "r", encoding="utf-8") as f:
        
        #Read lines
        for ln, raw in enumerate(f, 1):
            #Split to items
            s = raw.strip()
            parts = s.split()
            
            #Not enough columns
            if len(parts) < col - 1:
                continue
            
            #Append index to the list
            try:
                ids.append(int(float(parts[col])))
            
            #Impossible, skip
            except ValueError:
                continue
    return ids


def write_output_coordinates(out_path: str, centroid_ids: List[int], clusters: Dict[int, List[int]], points: Dict[int, Coord]):
    #Write coordinates of centroids and connected clients
    with open(out_path, "w", encoding="utf-8") as out:
        #Process centroids
        for cid in centroid_ids:
            #Get its coordinates
            centr_coord = points.get(cid)
            
            #Coordinates do not exist in the list
            if centr_coord is None:
                    continue
            
            #Write point
            else:
                    #print(cid, ':')
                    x, y, z = centr_coord 
                    out.write(f"{cid}\t{x}\t{y}\t{z}\n")
            
            #Get all connected points
            pt_ids = clusters.get(cid)
            
            #No connected points
            if not pt_ids:
                continue
            
            #Browse all connected points
            for pid in pt_ids:
                #print('    ', pid)
                #Get point coordinates
                coord = points.get(pid)
                
                #Coordinates do not exist in the list
                if coord is None:
                    continue
                
                #Write point
                else:
                    x, y, z = coord
                    out.write(f"{pid}\t{x}\t{y}\t{z}\n")



def main(argv: List[str] | None = None):
  
    parser = argparse.ArgumentParser(
        description="Collect coordinates of points connected to selected centroids."
    )
    
    parser.add_argument(
        "--points",
        required=True,
        help="Path to points file (columns: x y z).",
    )
    
    parser.add_argument("--clusters", required=True, 
        help="Path to clusters file (centroid_id: followed by point ids).",
    )
    
    parser.add_argument( "--centroids", required=True,
        help="Path to centroid ids file (one id per line).",
    )
    
    parser.add_argument("--column", type=int, required=True,
        help="Column with point ids.",
    )
    
    parser.add_argument(
        "--out",
        required=True,
        help="Path to write the grouped result file.",
    )
        
    print('\n>> Exporting clusters and their clients: \n')

    #Get parser arguments
    args = parser.parse_args(argv)

    print('Reading points')
    points = load_points(args.points)
    
    print('Reading clusters')
    clusters = load_clusters(args.clusters)
    
    print('Reading centroids')
    centroid_ids = load_centroid_ids(args.centroids, args.column)
    

    #if not centroid_ids:
    #    print("Warning: no centroid ids found.", file=sys.stderr)

    print('Writing results')
    write_output_coordinates(args.out, centroid_ids, clusters, points)
    
    return 0


if __name__ == "__main__":
    main()