#!/usr/bin/env python3

# from modules import *
import psycopg2

#####################################################
##  Database Connect
#####################################################

'''
Connects to the database using the connection string
'''
def openConnection():
# connection parameters - ENTER YOUR LOGIN AND PASSWORD HERE
    userid = "y18s2c9120_unikey"
    passwd = "password"
    myHost = "soit-db-pro-2.ucc.usyd.edu.au"

    # Create a connection to the database
    conn = None
    try:
        # Parses the config file and connects using the connect string
        conn = psycopg2.connect(database=userid,
                                    user=userid,
                                    password=passwd,
                                    host=myHost)
    except psycopg2.Error as sqle:
        print("psycopg2.Error : " + sqle.pgerror)
    
    # return the connection to use
    return conn

'''
Q1: List all the user associated issues in the database for a given user 
See assignment description for how to load user associated issues based on the user id (user_id)
'''
# without stored procedure
def findUserIssues(user_id):

    try:
        curs = conn.cursor()
        # execute the query  #/
        # SQL Query
        query = """SELECT * 
                    FROM A3_ISSUE 
                    WHERE (RESOLVER=%s OR CREATOR=%s OR VERIFIER=%s) 
                    ORDER BY TITLE ASC"""
        # input parameters
        params = (user_id, user_id, user_id,)
        # connect to database & fetch the query result
        curs.execute(query, params)
        # Print out the input user id
        print("USER ID: " + str(user_id) + ", " + "TrackIssue regard to this user is...\n")
        #  loop through the resultset to output tracker issues
        nr = 0
        row = curs.fetchone()
        while row is not None:

                nr+=1
                print("Issue ID: " + str(row[5]) + "\nTitle: " + str(row[0]) + "\nCreator: " + str(row[1]) + "\nResolver: " + str(row[2]) + "\nVerifier: " + str(row[3]) + "\nDescription " + str(row[4]) + "\n")
                row = curs.fetchone()

        if ( nr == 0 ):
            print("No entries found.")

        # clean up #/
        curs.close()

    except psycopg2.Error as sqle:       

        print("psycopg2.Error : " + sqle.pgerror)

# with stored procedure
def findUserIssues_storedProcesure(user_id):

    try:
        curs = conn.cursor()
        # execute the query  #/
        # SQL Query
        # connect to database & fetch the query result
        curs.callproc("LIST_A3ISSUE_BYUSERID", [user_id])
        # Print out the input user id
        print("USER ID: " + str(user_id) + ", " + "TrackIssue regard to this user is...\n")
        #  loop through the resultset to output tracker issues
        nr = 0
        row = curs.fetchone()
        while row is not None:

                nr+=1
                print("Issue ID: " + str(row[5]) + "\nTitle: " + str(row[0]) + "\nCreator: " + str(row[2]) + "\nResolver: " + str(row[3]) + "\nVerifier: " + str(row[4]) + "\nDescription " + str(row[1]) + "\n")
                row = curs.fetchone()

        if ( nr == 0 ):
            print("No entries found.")

        # clean up #/
        curs.close()

    except psycopg2.Error as sqle:       

        print("psycopg2.Error : " + sqle.pgerror)

'''
Q2: Find the associated issues for the user with the given userId (user_id) based on the searchString provided as the parameter, and based on the assignment description
'''
def findIssueBasedOnExpressionSearchedOnTitleAndDescription(searchString, user_id):

    try:
        curs = conn.cursor()
        # execute the query  #/
        query = """SELECT * 
                   FROM A3_ISSUE 
                   WHERE ((RESOLVER=%s OR CREATOR=%s OR VERIFIER=%s) AND (DESCRIPTION LIKE %s OR TITLE LIKE %s)) 
                   ORDER BY TITLE ASC"""
        # input parameters
        params = (user_id, user_id, user_id, "%" + searchString + "%", "%" + searchString + "%",)
        # connect to database & fetch the query result
        curs.execute(query, params)
        # Print out the input user id
        print("USER ID: " + str(user_id) + ", " + "Search keywords: " + searchString + "\n")
        #  loop through the resultset to output tracker issues
        nr = 0
        row = curs.fetchone()
        while row is not None:

                nr+=1
                print("Issue ID: " + str(row[5]) + "\nTitle: " + str(row[0]) + "\nCreator: " + str(row[1]) + "\nResolver: " + str(row[2]) + "\nVerifier: " + str(row[3]) + "\nDescription " + str(row[4]) + "\n")
                row = curs.fetchone()

        if ( nr == 0 ):
            print("No entries found.")

        # clean up #/
        curs.close()

    except psycopg2.Error as sqle:       

        print("psycopg2.Error : " + sqle.pgerror)

'''
Q3: Issue (new_issue, get all, get details, Add the details for a new issue to the database - details for new issue provided as parameters
'''
def addIssue(title, creator, resolver, verifier, description):
    
    try:
        curs = conn.cursor()
        # execute the query  #/
        query = """INSERT INTO A3_ISSUE (TITLE,DESCRIPTION,CREATOR,RESOLVER,VERIFIER) VALUES (%s,%s,%s,%s,%s)"""
        # input parameters
        params = (title, description, creator, resolver, verifier,)
        # connect to database & insert the new tracking record
        curs.execute(query, params)
        # Print out the input user id
        print("Insert new tracking issue... \n")
        conn.commit()
        curs.close()

    except psycopg2.Error as sqle:  
        # return False if succeed in SQL update     
        output = False
        print("psycopg2.Error : " + sqle.pgerror)
    
    else:
        # return True if succeed in SQL update
        output = True
    
    # TODO - add an issue
    # Insert a new issue to database
    # return False if adding was unsuccessful 
    # return Ture if adding was successful

    return output

'''
Q4: Update the details of an issue having the provided issue_id with the values provided as parameters
'''
def updateIssue(title, creator, resolver, verifier, description, issue_id):
     
    try:
        curs = conn.cursor()
        # execute the query  #/
        query = """
                UPDATE A3_ISSUE
                SET TITLE=%s, DESCRIPTION=%s, CREATOR=%s, RESOLVER=%s, VERIFIER=%s
                WHERE ISSUE_ID=%s
                """
        # input parameters
        params = (title, description, creator, resolver, verifier, issue_id,)
        # connect to database & insert the new tracking record
        curs.execute(query, params)
        # Print out the input user id
        print("Update new tracking issue... \n")
        conn.commit()
        curs.close()

    except psycopg2.Error as sqle:  
        # return False if succeed in SQL update   
        output = False
        print("psycopg2.Error : " + sqle.pgerror)   
    else:
        # return True if succeed in SQL update
        output = True

    return output

##
# Main program.
#/
'''
Function 1 - viewing issue list
''' 
# connect to database
conn = openConnection() 
## Function 1.1 without using STORED PROCEDURES
print("\n*** Question 1 PART 1 FindUserIssues -- Function without using Stored Procedures ***")

# Scenario 1: list user id = 1
user_id = 1 # define input values
findUserIssues(user_id) # apply function

# Scenario 2: list user id = 2
user_id = 2 # define input values
findUserIssues(user_id) # apply function

# Scenario 3: list user id = 3
user_id = 3 # define input values
findUserIssues(user_id) # apply function

## Function 1 Using STORED PROCEDURES
print("\n*** Question 1 PART 2 FindUserIssues -- Function using Stored Procedures ***")
# Scenario 1: list user id = 1
user_id = 1 # define input values
findUserIssues_storedProcesure(user_id) # apply function

# Scenario 2: list user id = 2
user_id = 2 # define input values
findUserIssues_storedProcesure(user_id) # apply function

# Scenario 3: list user id = 3
user_id = 3 # define input values
findUserIssues_storedProcesure(user_id) # apply function

'''
Function 2 - Search to Find Issues
'''
print("\n*** Question 2 findIssueBasedOnExpressionSearchedOnTitleAndDescription ***")
# Scenario 1: list keyword: zero
# define input values
searchString = 'zero'
user_id = 1
# apply function
findIssueBasedOnExpressionSearchedOnTitleAndDescription(searchString, user_id)
# Scenario 2: list no keyword
# define input values
searchString = ''
user_id = 1
# apply function
findIssueBasedOnExpressionSearchedOnTitleAndDescription(searchString, user_id)
# Scenario 3: list keyword: harsh sign
# define input values
searchString = '#'
user_id = 1
# apply function
findIssueBasedOnExpressionSearchedOnTitleAndDescription(searchString, user_id)

'''
Function 3 - Add new record in database
'''
print("\n*** Question 3 addIssue ***")
# Scenario 1: Succeed in Add-in new issue for existing user-id
# define input values
title = 'test'
creator = 1
resolver = 1
verifier = 1
description = 'test'
# apply function
output = addIssue(title, creator, resolver, verifier, description)
print("Scenario 1 Add-in new issue for existing user-id: " + str(output))

# Scenario 2: Add-in other Issue
# define input values
title = 'test 2'
creator = 2
resolver = 2
verifier = 1
description = 'test 2'
# apply function
output = addIssue(title, creator, resolver, verifier, description)
print("Scenario 2 Add-in result: " + str(output))

'''
Function 4 - Update record in database
'''
print("\n*** Question 4 updateIssue ***")
# Scenario 1: Succeed in Update exist user -id issue
# define input values
title = 'Division by zero'
creator = 1
resolver = 1
verifier = 1
issue_id = 2
description = 'Division by 0 doesn\'t yield error or infinity as would be expected. Instead it results in -1.'
# apply function
output = updateIssue(title, creator, resolver, verifier, description, issue_id)
print("Scenario 1 Update result: " + str(output))

# Scenario 2: Other update operation
title = 'Update test'
creator = 2
resolver = 2
verifier = 1
issue_id = 3
description = '#'
# apply function
output = updateIssue(title, creator, resolver, verifier, description, issue_id)
print("Scenario 2 Update result: " + str(output))