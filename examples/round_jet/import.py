from neo4j import GraphDatabase

'''
Steps to take atomization data from a merge_split.csv file created by ASSET
1. Start Neo4j and create a new database and local dbms
2. Move merge_split.csv file into the 'import' directory in your new dbms
3. Activate APOC in your DBMS
4. Start the dbms
5. Rename the relevant parameters below
6. Run the script
'''

# Connect to your Neo4j database !!! RENAME DEPENDING ON YOUR DBMS !!!
uri = "neo4j://localhost:7687"
username = "neo4j"
password = "12345678"  # replace with your password

# Create driver instance
driver = GraphDatabase.driver(uri, auth=(username, password))

# Function to execute a Cypher query
def execute_query(session, query):
    session.run(query)
    print(f"Query executed: {query[:50]}...")

# Step 1: Load CSV and create nodes 
# !!! ADJUST PARAMETERS DEPENDING ON YOUR MERGE_SPLIT FILE !!!
def load_csv_and_create_nodes(session):
    query = """
    CALL {
LOAD CSV WITH HEADERS FROM "file:///merge_split_clean.csv" AS csvline
    CREATE (n:droplet { 
        id: toInteger(csvline.NewID),
        Event: csvline.EventType,
        OldIDs: [x IN split(csvline.OldIDs, ";") WHERE x <> '' | toInteger(x)],
        NewID: toInteger(csvline.NewID),
        Volume: toFloat(csvline.NewVol),
        Event_Time: toFloat(csvline.Time),
        X: toFloat(csvline.X),
        Y: toFloat(csvline.Y),
        Z: toFloat(csvline.Z),
        U: toFloat(csvline.U),
        V: toFloat(csvline.V),
        W: toFloat(csvline.W),
        L1: toFloat(csvline.L1),
        L2: toFloat(csvline.L2),
        L3: toFloat(csvline.L3)
    })
    } IN TRANSACTIONS OF 10000 ROWS;
    """
    execute_query(session, query)

# Step 2: Create Merge Relations
def create_merge_relations(session):
    query1 = """
    MATCH (m:droplet {Event: "Merge"})
    WITH m, [id IN m.OldIDs WHERE id IS NOT NULL] AS oldIDs
    WHERE NONE(id IN oldIDs WHERE id = 1)
    UNWIND oldIDs AS oldID
    MATCH (d:droplet {NewID: oldID})
    CREATE (d)-[:Merge]->(m)
    """
    execute_query(session, query1)

    query2 = """
    MATCH (m:droplet {Event: "Merge"})
    WITH m, [id IN m.OldIDs WHERE id IS NOT NULL] AS oldIDs
    WHERE ANY(id IN oldIDs WHERE id = 1)
    UNWIND oldIDs AS oldID
    MATCH (d:droplet {NewID: oldID}),(c:droplet {Event:"None"})
    Where d.NewID <> 1
    CREATE (d)-[:Merge]->(c)
    """
    execute_query(session, query2)

# Step 3: Delete Merge Nodes
def delete_merge_nodes(session):
    query = """
    MATCH (m:droplet {Event: "Merge"})
    WHERE m.NewID = 1
    delete m
    """
    execute_query(session, query)

# Step 4: Create Split Relations
def create_split_relations(session):
    query = """
    MATCH (n:droplet {Event: "Split"}), (d:droplet)
    WHERE ANY(id IN n.OldIDs WHERE id = d.NewID)
    CREATE (d)-[:Split]->(n)
    """
    execute_query(session, query)

# Step 5: Delete "Out of Domain" Nodes
def delete_out_of_domain_nodes(session):
    query = """
    MATCH (n:droplet)
    WHERE n.Event = "Out of Domain"
    DETACH DELETE n
    """
    execute_query(session, query)

# Function to rename nodes
def rename_labels(session, label_match, label_new, condition):
    query = f"""
    MATCH (n:{label_match})
    WHERE {condition}
    WITH collect(n) as p
    CALL apoc.refactor.rename.label("{label_match}","{label_new}",p)
    YIELD committedOperations
    RETURN committedOperations
    """
    execute_query(session, query)

# Function to rename remaining droplet nodes
def rename_remaining_nodes(session, current_label, next_label):
    query = f"""
    MATCH (n:droplet),(d:`{current_label}`)
    WHERE ANY(id IN n.OldIDs WHERE id = d.NewID)
    WITH collect(n) as p
    CALL apoc.refactor.rename.label("droplet","{next_label}",p)
    YIELD committedOperations
    RETURN committedOperations
    """
    execute_query(session, query)

# Main
with driver.session() as session:
    load_csv_and_create_nodes(session)
    create_merge_relations(session)
    delete_merge_nodes(session)
    create_split_relations(session)
    delete_out_of_domain_nodes(session)
    
    # Renaming labels
    rename_labels(session, "droplet", "core", 'n.Event = "None"')
    rename_labels(session, "droplet", "1", "ANY(id IN n.OldIDs WHERE id = 1)")

    # Rename remaining "droplet" nodes
    current_label = 1
    while True:
        droplet_count = session.run("MATCH (n:droplet) RETURN count(n) AS count").single()[0]
        if droplet_count == 0:
            print("No more 'droplet' nodes to rename. Process completed.")
            break

        next_label = current_label + 1
        rename_remaining_nodes(session, current_label, next_label)
        current_label += 1

# Close the driver connection
driver.close()
