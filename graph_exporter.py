import pandas as pd

def generate_cypher_queries(df: pd.DataFrame, output_file: str = "knowledge_graph.cypher"):
    """
    Translates the final dataframe into Neo4j Cypher queries for ingestion.
    Creates Chemical nodes, Cluster nodes, and relationships.
    """
    queries = []
    
    # Ensure Constraints
    queries.append("CREATE CONSTRAINT IF NOT EXISTS FOR (c:Chemical) REQUIRE c.smiles IS UNIQUE;")
    queries.append("CREATE CONSTRAINT IF NOT EXISTS FOR (cl:Cluster) REQUIRE cl.id IS UNIQUE;")
    queries.append("CREATE CONSTRAINT IF NOT EXISTS FOR (t:Target) REQUIRE t.name IS UNIQUE;")
    
    clusters_added = set()
    targets_added = set()
    
    for _, row in df.iterrows():
        smiles = row.get('Standardized SMILES', '')
        if not smiles or pd.isna(smiles):
            continue
            
        cluster_id = row.get('Cluster ID', 'Unclustered')
        decision = row.get('Decision', 'Unknown')
        domain = row.get('Domain', 'Unknown')
        
        # 1. Create Cluster Node
        if cluster_id not in clusters_added:
            justification = row.get('Read_Across_Justification', '')
            # Escape quotes
            justification = str(justification).replace('"', "'")
            c_query = f'MERGE (cl:Cluster {{id: "{cluster_id}"}}) SET cl.decision = "{decision}", cl.justification = "{justification}"'
            queries.append(c_query)
            clusters_added.add(cluster_id)
            
        # 2. Create Chemical Node
        chem_query = f'MERGE (c:Chemical {{smiles: "{smiles}"}}) SET c.domain = "{domain}"'
        queries.append(chem_query)
        
        # 3. Create BELONGS_TO relationship
        rel_query = f'MATCH (c:Chemical {{smiles: "{smiles}"}}), (cl:Cluster {{id: "{cluster_id}"}}) MERGE (c)-[:BELONGS_TO]->(cl)'
        queries.append(rel_query)
        
        # 4. Create Biological Target Nodes and SHARES_TARGET relationships
        targets_str = row.get('ChEMBL_Targets', '')
        if pd.notna(targets_str) and targets_str != "None":
            targets = [t.strip() for t in str(targets_str).split(',')]
            for t in targets:
                t_clean = t.replace('"', "'")
                if t_clean not in targets_added:
                    queries.append(f'MERGE (t:Target {{name: "{t_clean}"}})')
                    targets_added.add(t_clean)
                queries.append(f'MATCH (c:Chemical {{smiles: "{smiles}"}}), (t:Target {{name: "{t_clean}"}}) MERGE (c)-[:HAS_ACTIVITY_AGAINST]->(t)')
                
    with open(output_file, 'w', encoding='utf-8') as f:
        for q in queries:
            f.write(q + ';\n')
            
    return output_file
