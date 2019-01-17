# Importing the frameworks

from modules import *
from flask import *
import database

user_details = {}                   # User details kept for us
session = {}
page = {}

# Initialise the application
app = Flask(__name__)
app.secret_key = 'aab12124d346928d14710610f'


#####################################################
##  INDEX
#####################################################

@app.route('/')
def index():
    # Check if the user is logged in
    if('logged_in' not in session or not session['logged_in']):
        return redirect(url_for('login'))
    page['title'] = 'IssueTracker'
    
    return redirect(url_for('list_issue'))

    #return render_template('index.html',
    #    session=session,
    #    page=page,
    #    user=user_details)

#####################################################
##  LOGIN
#####################################################

@app.route('/login', methods=['POST', 'GET'])
def login():
    # Check if they are submitting details, or they are just logging in
    if(request.method == 'POST'):
        # submitting details
        login_return_data = check_login(request.form['id'])

        # If it's null, saying they have incorrect details
        if login_return_data is None:
            page['bar'] = False
            flash("Incorrect id, please try again")
            return redirect(url_for('login'))

        # If there was no error, log them in
        page['bar'] = True
        flash('You have been logged in successfully')
        session['logged_in'] = True

        # Store the user details for us to use throughout
        global user_details
        user_details = login_return_data
        return redirect(url_for('index'))

    elif(request.method == 'GET'):
        return(render_template('login.html', page=page))

#####################################################
##  LOGOUT
#####################################################

@app.route('/logout')
def logout():
    session['logged_in'] = False
    page['bar'] = True
    flash('You have been logged out')
    return redirect(url_for('index'))

#####################################################
##  LIST ISSUE
#####################################################

@app.route('/issue', methods=['POST', 'GET'])
def list_issue():
    if( 'logged_in' not in session or not session['logged_in']):
        return redirect(url_for('login'))
    # The user is just viewing the page
    if (request.method == 'GET'):
        # First check if specific event
        issue_list = database.findUserIssues(user_details['user_id'])
        if(issue_list is None):
            issue_list = []
            flash("There are no issues in our system for user " + user_details['user_id'])
            page['bar'] = False
        return render_template('issue_list.html', issue=issue_list, session=session, page=page)

    # Try to get from the database
    elif(request.method == 'POST'):
        issue_list_find = database.findIssueBasedOnExpressionSearchedOnTitleAndDescription(request.form['search'], user_details['user_id'])
        if(issue_list_find is None):
            issue_list_find = []
            flash("Issue \'{}\' does not exist for user ".format(request.form['search']) + user_details['user_id'])
            page['bar'] = False
        return render_template('issue_list.html', issue=issue_list_find, session=session, page=page)

#####################################################
##  Add Issue
#####################################################

@app.route('/new-issue' , methods=['GET', 'POST'])
def new_issue():
    if( 'logged_in' not in session or not session['logged_in']):
        return redirect(url_for('login'))

    # If we're just looking at the 'new issue' page
    if(request.method == 'GET'):
        times = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
        return render_template('new_issue.html', user=user_details, times=times, session=session, page=page)
    # If we're adding a new issue
    success = database.addIssue(request.form['title'],
                                 request.form['creator'],
                                 request.form['resolver'],
                                 request.form['verifier'],
                                 request.form['description'])
    if(success == True):
        page['bar'] = True
        flash("Issue Added!")
        return(redirect(url_for('index')))
    else:
        page['bar'] = False
        flash("There was an error adding new issue.")
        return(redirect(url_for('new_issue')))

#####################################################
## UPDATE ISSUE
#####################################################
@app.route('/update_issue/', methods=['GET', 'POST'])
def update_issue():
    if( 'logged_in' not in session or not session['logged_in']):
        return redirect(url_for('login'))

    # Check the details of the issue
    issue_id = request.args.get('issue_id')

    #if not issue_id:
    #    page['bar'] = False
    #    flash("Error, no issue was given. URL requires \'?issue_id=<id>\'")
    #    return(redirect(url_for('index')))

    issue_results = get_issue(issue_id,user_details['user_id'])
    #print(issue_results)

    if issue_results is None:
        issue_results = []

    # If we're just looking at the 'update issue' page
    if(request.method == 'GET'):
        times = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
        return render_template('update_issue.html', issueInfo = issue_results, user=user_details, times=times, session=session, page=page)
    # If we're updating an issue
    success = database.updateIssue(request.form['title'],
                                 request.form['creator'],
                                 request.form['resolver'],
                                 request.form['verifier'],
                                 request.form['description'],
                                 request.form['issue_id'])
    if(success == True):
        page['bar'] = True
        flash("Issue Updated!")
        return(redirect(url_for('index')))
    else:
        page['bar'] = False
        flash("There was an error updating the issue.")
        return(redirect(url_for('update_issue')))
    
def get_issue(issue_id, user_id):
    print('routes.getIssue')
    for issue in database.findUserIssues(user_id):
        #print(issue['issue_id'])
        #print(issue_id)
        if issue['issue_id'] == issue_id:
            return [issue]
    return None	

def check_login(user_id):
    print('routes.check_login')
    user_data = [user_id, '', '']

    tuples = {
            'user_id': user_data[0],
            'first_name': user_data[1],
            'family_name': user_data[2],
        }
    return tuples
