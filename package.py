from pycman.utils.common import datetime_format

package = {
    'name': 'biokit',
    'version': '0.0.1',
    'author': 'tanghongzhen',
    'email': 'tanghongzhen34@gmail.com',
    'scripts': {
        'default': 'echo hello!',
        'tests': f'pytest tests -n=auto --html=test-reports/test-report-{datetime_format()}.html --self-contained-html'
    }
}