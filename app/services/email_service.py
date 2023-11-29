import os
import smtplib
import jinja2
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.header import Header
from email.utils import formataddr

class SingleEmailHeader(object):
    def __init__(self, recipients, email_params, template=''):
        self.recipients = recipients
        self.server = os.environ.get("EMAIL_SERVER")
        self.from_email = os.environ.get("EMAIL_FROM_EMAIL")
        self.from_name = os.environ.get("EMAIL_FROM_NAME")
        ## Construct email o
        self.msg = MIMEMultipart('alternative')
        self.msg['Subject'] = email_params['subject']
        self.msg['From'] = formataddr((str(Header(self.from_name, 'utf-8')), self.from_email))
        self.msg['To'] = self.recipients
        ## Render email HTML content
        # self.html = self.render(os.path.join(os.path.dirname(__file__), 'templates', template), email_params)
        # self.msg.attach(MIMEText(self.html, 'html'))
        message_content_text = email_params['message']
        message_content = MIMEText(message_content_text, "plain")
        self.msg.attach(message_content)

    def render(self, tpl_path, email_params):
        path, filename = os.path.split(tpl_path)
        return jinja2.Environment(loader=jinja2.FileSystemLoader(path or './')).get_template(filename).render(email_params)

    def sendmail(self):
        smtp_server = smtplib.SMTP(self.server, port=25)
        # The TO and CC header fields are populated by the header construction, and any additional recipient addresses are effectively BCC
        smtp_server.sendmail(self.from_email, self.recipients, self.msg.as_string())
        smtp_server.quit()

class EmailService:
    def __init__(self) -> None:
        pass

    def send_email(self, recipients=None, email_subject='', message_body='', template_file=None):
        if isinstance(recipients, list):
            recipients = recipients[0]
        email_params = {
            "subject": email_subject,
            "message": message_body,
        }
        try:
            email = SingleEmailHeader(recipients, email_params, template=template_file)
        except Exception as e:
            print("Error creating header:")
            print(e)
        try:
            email.sendmail()
        except Exception as e:
            print("Error sending email:")
            print(e)